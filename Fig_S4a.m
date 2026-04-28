R = [100:50:250 300:100:800 1000];
R_all = []; 
theta_all = []; 
period_all = []; 
theta_interp = 0.005:0.005:10;
for i = 1:9
    E_MD = readmatrix("energy_R"+num2str(R(i))+".dat");
    theta = E_MD(:,1);
    energy = E_MD(:,2) - E_MD(1,2);

    f_energy = spline(theta,energy);
    energy_interp = ppval(f_energy,theta_interp);

    [energy_min,loc_min] = findpeaks(-energy_interp,'MinPeakDistance',10,'MinPeakWidth',4);
    theta_min = theta_interp(loc_min);

    period_temp = diff(theta_min);
    period_avg = mean(period_temp);
    index = period_temp < period_avg*1.2 & period_temp > period_avg*0.8;

    period_all = [period_all period_temp(index)];
    theta_all = [theta_all theta_min(index)];
    R_all = [R_all, R(i)*ones(1,sum(index))];
end
%% R800
E_MD = readmatrix("energy_R"+num2str(800)+".dat");
theta = E_MD(:,1); energy = E_MD(:,2) - E_MD(1,2);
f_energy_2 = spline(theta(400:end),energy(400:end));
theta_interp_1 = theta(1:400)';
theta_interp_2 = 0.4:0.005:10;
theta_interp = [theta_interp_1 theta_interp_2];
energy_interp_1 = energy(1:400)';
energy_interp_2 = ppval(f_energy_2,theta_interp_2);
[energy_min_1,loc_min_1] = findpeaks(-energy_interp_1);
[energy_min_2,loc_min_2] = findpeaks(-energy_interp_2,'MinPeakDistance',5,'MinPeakWidth',3);
energy_min = [energy_min_1 energy_min_2];
loc_min = [loc_min_1 loc_min_2+400];

theta_min = theta_interp(loc_min);
period_temp = diff(theta_min);
period_avg = mean(period_temp);
index = period_temp < period_avg*1.2 & period_temp > period_avg*0.8;

period_all = [period_all period_temp(index)];
theta_all = [theta_all theta_min(index)];
R_all = [R_all, 800*ones(1,sum(index))];

%% R1000
E_MD = readmatrix("energy_R"+num2str(1000)+".dat");
theta = E_MD(:,1); energy = E_MD(:,2) - E_MD(1,2);
f_energy_2 = spline(theta(400:end),energy(400:end));
theta_interp_1 = theta(1:400)';
theta_interp_2 = 0.4:0.005:10;
theta_interp = [theta_interp_1 theta_interp_2];
energy_interp_1 = energy(1:400)';
energy_interp_2 = ppval(f_energy_2,theta_interp_2);
[energy_min_1,loc_min_1] = findpeaks(-energy_interp_1);
[energy_min_2,loc_min_2] = findpeaks(-energy_interp_2,'MinPeakDistance',5,'MinPeakWidth',3);
energy_min = [energy_min_1 energy_min_2];
loc_min = [loc_min_1 loc_min_2+400];

theta_min = theta_interp(loc_min);
period_temp = diff(theta_min);
period_avg = mean(period_temp);
index = period_temp < period_avg*1.2 & period_temp > period_avg*0.8;

period_all = [period_all period_temp(index)];
theta_all = [theta_all theta_min(index)];
R_all = [R_all, 1000*ones(1,sum(index))];


%% Fit A/x  0.9991
% A = 56.9749 (56.9026, 57.0471)
figure; hold on
scatter(R_all,period_all,'k.');
x = 100:1000;
plot(x,56.9749./x,'k-')
xlabel("Flake size",'Interpreter','latex')
ylabel("Period ($^\circ$)",'Interpreter','latex')