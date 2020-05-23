import numpy as np
from numpy import linalg as LA
from PIL import Image
import spacecraft_params as sc
import const as cn

cos    = np.cos
sin    = np.sin
sqrt   = np.sqrt
tanh   = np.tanh
cross  = np.cross
matmul = np.matmul
norm   = LA.norm


__all__ = ["mee2rv", "inertial2radial", "thrust_angle", "drawEarth"]

################################################################################
def mee2rv(p,f,g,h,k,L,mu):
    """ Convert modified equinoctial elements to cartesian coordinates
    Robyn Woollands 04/13/2020
    Parameters:
    p  -- semilatus rectum
    f  -- e*cos(w+Omega)
    g  -- e*sin(w+Omega)
    h  -- tan(i/2)*cos(Omega)
    k  -- tan(i/2)*sin(Omega)
    L  -- true longitude
    mu -- standard gravitational parameter
    Returns:
    r -- position
    v -- velocity
    """
    # Common terms
    alpha2 = h**2 - k**2
    tani2s = h**2 + k**2
    s2     = 1 + tani2s
    cosL   = cos(L)
    sinL   = sin(L)
    hk2    = 2*h*k
    sq     = sqrt(mu/p)
    # Radius
    radius = p/(1 + f*cosL + g*sinL)
    # Position and Velocity
    r = np.zeros((np.size(p),3))
    v = np.zeros((np.size(p),3))
    r[:,0] = radius*(cosL + alpha2*cosL + hk2*sinL)/s2
    r[:,1] = radius*(sinL - alpha2*sinL + hk2*cosL)/s2
    r[:,2] = 2*radius*(h*sinL - k*cosL)/s2
    v[:,0] = -sq*(sinL + alpha2*sinL - hk2*cosL + g - f*hk2 + alpha2*g)/s2
    v[:,1] = -sq*(-cosL + alpha2*cosL + hk2*sinL - f + g*hk2 + alpha2*f)/s2
    v[:,2] = 2*sq*(h*cosL + k*sinL + f*h + g*k)/s2
    return r, v
################################################################################

def inertial2radial(r,v):

    # Radial in x- and z- direction
    hvec = cross(r,v);
    h    = norm(hvec);
    xrdl = r/norm(r);
    zrdl = hvec/h;

    # Radial in y-direction
    yrdl = cross(zrdl,xrdl);

    return np.array([xrdl,yrdl,zrdl])

################################################################################
def thrust_angle(data,rho,eclipse):

    ind     = np.size(data,0)
    S       = np.zeros(ind)
    F       = np.zeros(ind)
    delta   = np.zeros(ind)
    zeta    = np.zeros(ind)
    Pa      = np.zeros(ind)
    u_inert = np.zeros((ind,3))
    u_lvlh  = np.zeros((ind,3))

    for i in range(ind):

        p = data[i,1]
        f = data[i,2]
        g = data[i,3]
        h = data[i,4]
        k = data[i,5]
        L = data[i,6]
        m = data[i,7]
        plam = data[i,8]
        flam = data[i,9]
        glam = data[i,10]
        hlam = data[i,11]
        klam = data[i,12]
        Llam = data[i,13]
        mlam = data[i,14]

        [[r],[v]] = mee2rv(p,f,g,h,k,L,cn.mu)

        # Common terms
        SinL = sin(L)
        CosL = cos(L)
        q    = 1+f*CosL+g*SinL
        s    = 1+h**2+k**2
        C1   = sqrt(p)
        C2   = 1/q
        C3   = h*SinL-k*CosL

        B = np.array([[0,2*p*C2*C1,0],
        [C1*SinL,C1*C2*((q+1)*CosL+f),-C1*(g/q)*C3],
        [-C1*CosL,C1*C2*((q+1)*SinL+g),C1*(f/q)*C3],
        [0,0,C1*s*CosL*C2/2],
        [0,0,C1*s*SinL*C2/2],
        [0,0,C1*C2*C3],
        [0,0,0]])

        lam = np.array([plam,flam,glam,hlam,klam,Llam,mlam])
        BTL = matmul(B.T,lam)

        S[i]     = sc.c*sc.si2can*norm(BTL)/m+mlam-1
        delta[i] = 0.5*(1+tanh(S[i]/rho))

        # Eclipse model
        if (eclipse):
            if (r[0] < 0): # Assume Sun is located along positive x-axis
                gamma   = norm(r[1:2]) - cn.Req/cn.DU
                zeta[i] = 0.5*(1.0+tanh(gamma/rho))
                Pa[i]   = zeta[i]*sc.P
                Thr     = Pa[i]*sc.A*sc.eta/sc.Isp/cn.g0
                F[i]    = Thr*delta[i]
                u_inert[i,:] = -BTL/norm(BTL)*delta[i]*zeta[i]
            else:
                Pa[i] = sc.P
                Thr   = Pa[i]*sc.A*sc.eta/sc.Isp/cn.g0
                F[i]  = Thr*delta[i]
                u_inert[i,:] = -BTL/norm(BTL)*delta[i]
        else:
            Pa[i] = sc.P
            Thr   = sc.P*sc.A*sc.eta/sc.Isp/cn.g0
            F[i]  = Thr*delta[i]
            u_inert[i,:] = -BTL/norm(BTL)*delta[i]

        M = inertial2radial(r,v)
        u_lvlh[i,:] = matmul(M,u_inert[i,:].T).T

    return u_inert, u_lvlh, S, F, Pa, delta, zeta

################################################################################
def drawEarth(Radius):

    # Create a sphere with earths surface texture

    # Load texture
    #response = requests.get('http://www.johnstonsarchive.net/spaceart/cmaps/earthmap.jpg')

    #img = Image.open(StringIO(response.content))
    img = Image.open('blue_marble.jpg')

    # Rescale RGB values
    img = np.array(img.resize([int(d/4) for d in img.size]))/256

    # Image coordinates
    lons = np.linspace(-180, 180, img.shape[1]) * np.pi/180
    lats = np.linspace(-90, 90, img.shape[0])[::-1] * np.pi/180

    x = Radius*np.outer(np.cos(lons), np.cos(lats)).T
    y = Radius*np.outer(np.sin(lons), np.cos(lats)).T
    z = Radius*np.outer(np.ones(np.size(lons)), np.sin(lats)).T

    # Alternatively, create a simple sphere object (faster)
    # pi = np.pi
    # phi, theta = np.mgrid[0.0:pi:100j, 0.0:2.0*pi:100j]
    # x = Radius*np.sin(phi)*np.cos(theta)
    # y = Radius*np.sin(phi)*np.sin(theta)
    # z = Radius*np.cos(phi)

    return x, y, z, img
