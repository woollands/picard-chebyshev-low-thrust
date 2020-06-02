% NOTE: The following code (17 lines) was obtained from MathWorks
% online file exchange (Ryan Gray).
load('topo.mat','topo','topomap1');
colormap(topomap1);
% Create the surface.
[x,y,z] = sphere(50);
props.AmbientStrength           = 0.1;
props.DiffuseStrength           = 1;
props.SpecularColorReflectance  = .5;
props.SpecularExponent          = 20;
props.SpecularStrength          = 1;
props.FaceColor                 = 'texture';
props.EdgeColor                 = 'none';
props.FaceLighting              = 'phong';
props.Cdata                     = topo;
figure(1)
surface(x,y,z,props);
%     set(gca,'color','black')
axis equal