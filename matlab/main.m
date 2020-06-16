clear all; close all; clc;

const;
load_spice_kernels;

% Time
t0 = cspice_str2et('Jan 1 2020 00:00:00');
tf = cspice_str2et('Jan 1 2020 06:00:00');
time = [t0:1:tf];

% Initial Conditions (circular test orbit, no thrust)
[r0,v0] = elm2rv(1000+Req,0,0,0,0,0,0,mu,0);

% Compute analytical two-body solutions
for i = 1:length(time)
    [Rout(i,:),Vout(i,:)] = FnG(0,time(i),r0,v0,mu);
end

figure(1)
plot3(Rout(:,1)./Req,Rout(:,2)./Req,Rout(:,3)./Req,'k-','Linewidth',2)
hold on; grid on;
xlabel('X (km)')
ylabel('Y (km)')
zlabel('Z (km)')

% Eclipse Model
for i = 1:length(time)
    [TF(i),PA(i),ang(i),alpha(i)] = eclipse_grazinggoat(time(i),Rout(i,:));
end

figure(2)
plot((time-time(1))./3600,PA,'bo-')
ylabel('Power (W)')
xlabel('Time (hours)')
title('Power (Grazing Goat Model)')

ind1 = find(TF == 0);
ind2 = find(TF > 0 & TF < 1);

figure(1)
plot3(Rout(ind1,1)./Req,Rout(ind1,2)./Req,Rout(ind1,3)./Req,'r.','MarkerSize',30)
plot3(Rout(ind2,1)./Req,Rout(ind2,2)./Req,Rout(ind2,3)./Req,'g.','MarkerSize',30)

figure(1)
plotEarth;
legend('LEO','Umbra','Penumbra')

figure(3)
plot((time-time(1))./3600,ang.*180/pi)
hold on
plot((time(ind1)-time(1))./3600,ang(ind1).*180/pi,'ro')
title('Angle w.r.t. Penumbral Cone Axis')
xlabel('Time (hours)')
ylabel('Angle (deg)')

figure(4)
plot(time-time(1),TF)

