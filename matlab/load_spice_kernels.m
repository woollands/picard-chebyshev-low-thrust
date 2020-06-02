
pth = '/Users/robynmw/Dropbox/CODE/mice';
addpath( [pth,'/lib'])      
addpath( [pth,'/src/mice']) 
addpath( [pth,'/kernels'])

cspice_furnsh( [pth,'/kernels/de438.bsp'] )
cspice_furnsh( '/Users/robynmw/Dropbox/CODE/mice/kernels/naif0012.tls' )
