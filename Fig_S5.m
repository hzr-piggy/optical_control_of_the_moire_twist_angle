R = [100:50:250 300:100:800 1000];
R_all = []; R_small = []; R_large = [];
theta_all = []; theta_small = []; theta_large = [];
barrier_all = []; barrier_small = []; barrier_large = [];
theta_interp = 0.005:0.005:10;
for i = 1:9
    E_MD = readmatrix("energy_R"+num2str(R(i))+".dat");
    theta = E_MD(:,1);
    energy = E_MD(:,2) - E_MD(1,2);

    f_energy = spline(theta,energy);
    energy_interp = ppval(f_energy,theta_interp);

    [energy_min,loc_min] = findpeaks(-energy_interp,'MinPeakDistance',10,'MinPeakWidth',4);
    [energy_max,loc_max] = findpeaks(energy_interp,'MinPeakDistance',10,'MinPeakWidth',4);
    
    barrier = energy_max(1)+energy_min(1); loc_barrier = floor((loc_max(1)+loc_min(1))/2);
    shape_min = size(loc_min); shape_max = size(loc_max);
    n_max = 2; n_min = 2;
    while n_min <= shape_min(2) & n_max <= shape_max(2)
        if loc_min(n_min) > loc_max(n_max) & n_max < shape_max(2)
            n_max = n_max + 1;
        end
        if n_min < shape_min(2) & loc_min(n_min+1) < loc_max(n_max)
            n_min = n_min + 1;
        end
        barrier = [barrier energy_max(n_max)+energy_min(n_min)];
        loc_barrier = [loc_barrier floor((loc_max(n_max)+loc_min(n_min))/2)];
        n_min = n_min+1; n_max = n_max+1;
    end
    shape = size(barrier);
    theta_barrier = theta_interp(loc_barrier);
    R_all = [R_all, R(i)*ones(shape)];
    theta_all = [theta_all, theta_barrier];
    barrier_all = [barrier_all barrier];

    small_angle = theta_barrier<74.22/sqrt(R(i))-1;
    R_small = [R_small, R(i)*ones(1,sum(small_angle))];
    theta_small = [theta_small theta_barrier(small_angle)];
    barrier_small = [barrier_small barrier(small_angle)];
    % % 
    theta_large_temp = theta_barrier(~small_angle);
    barrier_large_temp = barrier(~small_angle);
    [barrier_large_max,loc_large] = findpeaks(barrier_large_temp);
    barrier_large = [barrier_large barrier_large_max];
    theta_large = [theta_large theta_large_temp(loc_large)];
    shape_large = size(loc_large);
    R_large = [R_large R(i)*ones(shape_large)];
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
[energy_max_1,loc_max_1] = findpeaks(energy_interp_1);
[energy_min_2,loc_min_2] = findpeaks(-energy_interp_2,'MinPeakDistance',5,'MinPeakWidth',3);
[energy_max_2,loc_max_2] = findpeaks(energy_interp_2,'MinPeakDistance',5,'MinPeakWidth',3);
energy_min = [energy_min_1 energy_min_2];
energy_max = [energy_max_1 energy_max_2];
loc_min = [loc_min_1 loc_min_2+400];
loc_max = [loc_max_1 loc_max_2+400];
barrier = energy_max(1)+energy_min(1); loc_barrier = floor((loc_max(1)+loc_min(1))/2);
shape_min = size(loc_min); shape_max = size(loc_max);
n_max = 2; n_min = 2;
while n_min <= shape_min(2) & n_max <= shape_max(2)
    if loc_min(n_min) > loc_max(n_max) & n_max < shape_max(2)
        n_max = n_max + 1;
    end
    if n_min < shape_min(2) & loc_min(n_min+1) < loc_max(n_max)
        n_min = n_min + 1;
    end
    barrier = [barrier energy_max(n_max)+energy_min(n_min)];
    loc_barrier = [loc_barrier floor((loc_max(n_max)+loc_min(n_min))/2)];
    n_min = n_min+1; n_max = n_max+1;
end
shape = size(barrier);
theta_barrier = theta_interp(loc_barrier);
R_all = [R_all, 800*ones(shape)];
theta_all = [theta_all, theta_barrier];
barrier_all = [barrier_all barrier];

small_angle = theta_barrier<74.22/sqrt(800)-1;
R_small = [R_small, 800*ones(1,sum(small_angle))];
theta_small = [theta_small theta_barrier(small_angle)];
barrier_small = [barrier_small barrier(small_angle)];
% % 
theta_large_temp = theta_barrier(~small_angle);
barrier_large_temp = barrier(~small_angle);
[barrier_large_max,loc_large] = findpeaks(barrier_large_temp);
barrier_large = [barrier_large barrier_large_max];
theta_large = [theta_large theta_large_temp(loc_large)];
shape_large = size(loc_large);
R_large = [R_large 800*ones(shape_large)];

%% R1000
E_MD = readmatrix("energy_R"+num2str(1000)+".dat");
theta = E_MD(:,1); energy = E_MD(:,2) - E_MD(1,2);
f_energy_2 = spline(theta(400:end),energy(400:end));
theta_interp_1 = theta(1:400)';
theta_interp_2 = 0.4:0.0025:10;
theta_interp = [theta_interp_1 theta_interp_2];
energy_interp_1 = energy(1:400)';
energy_interp_2 = ppval(f_energy_2,theta_interp_2);
[energy_min_1,loc_min_1] = findpeaks(-energy_interp_1);
[energy_max_1,loc_max_1] = findpeaks(energy_interp_1);
[energy_min_2,loc_min_2] = findpeaks(-energy_interp_2,'MinPeakDistance',5,'MinPeakWidth',3);
[energy_max_2,loc_max_2] = findpeaks(energy_interp_2,'MinPeakDistance',5,'MinPeakWidth',3);
energy_min = [energy_min_1 energy_min_2];
energy_max = [energy_max_1 energy_max_2];
loc_min = [loc_min_1 loc_min_2+400];
loc_max = [loc_max_1 loc_max_2+400];
barrier = energy_max(1)+energy_min(1); loc_barrier = floor((loc_max(1)+loc_min(1))/2);
shape_min = size(loc_min); shape_max = size(loc_max);
n_max = 2; n_min = 2;
while n_min <= shape_min(2) & n_max <= shape_max(2)
    if loc_min(n_min) > loc_max(n_max) & n_max < shape_max(2)
        n_max = n_max + 1;
    end
    if n_min < shape_min(2) & loc_min(n_min+1) < loc_max(n_max)
        n_min = n_min + 1;
    end
    barrier = [barrier energy_max(n_max)+energy_min(n_min)];
    loc_barrier = [loc_barrier floor((loc_max(n_max)+loc_min(n_min))/2)];
    n_min = n_min+1; n_max = n_max+1;
end
shape = size(barrier);
theta_barrier = theta_interp(loc_barrier);
R_all = [R_all, 1000*ones(shape)];
theta_all = [theta_all, theta_barrier];
barrier_all = [barrier_all barrier];

small_angle = theta_barrier<74.22/sqrt(1000)-1;
R_small = [R_small, 1000*ones(1,sum(small_angle))];
theta_small = [theta_small theta_barrier(small_angle)];
barrier_small = [barrier_small barrier(small_angle)];
% % 
theta_large_temp = theta_barrier(~small_angle);
barrier_large_temp = barrier(~small_angle);
[barrier_large_max,loc_large] = findpeaks(barrier_large_temp);
barrier_large = [barrier_large barrier_large_max];
theta_large = [theta_large theta_large_temp(loc_large)];
shape_large = size(loc_large);
R_large = [R_large 1000*ones(shape_large)];


figure; hold on
scatter3(R_all,theta_all,barrier_all,'k.')
scatter3(R_large,theta_large,barrier_large,20,'rs','filled')
zscale('log')

r_plot = 100:2:1100;
theta_plot = 0.02:0.02:11;
[R_plot,THETA_plot] = meshgrid(r_plot,theta_plot);
small_angle_plot = THETA_plot < 74.22./sqrt(R_plot);
large_angle_plot = THETA_plot >= 74.22./sqrt(R_plot);

theta_line = 0.02:0.02:11; shape_line = size(theta_line);
for i = 1:11
    small_angle_line = theta_line < 74.22/sqrt(R(i));
    large_angle_line = theta_line >= 74.22/sqrt(R(i));
    barrier_line = 1.4843e3*theta_line.^-3.*large_angle_line-0.0018*R(i)*sind(theta_line*3).*large_angle_line;
    plot3(R(i)*ones(shape_line),theta_line,barrier_line,'r--','LineWidth',1)
end

r_b = 100:1000;
theta_b = 74.22./sqrt(r_b);
E_b = 1.484e3*theta_b.^-3-0.0018*r_b.*sind(theta_b*3);
for n = 1:900
    fill3([r_b(n),r_b(n),r_b(n+1),r_b(n+1)], ...
        [theta_b(n),theta_b(n),theta_b(n+1),theta_b(n+1)], ...
        [1,E_b(n),E_b(n+1),1], 'r', 'FaceAlpha', 0.1, 'EdgeColor','none');
end

zlim([0.1,100])
ylim([2,10])
xlabel("Radius")
ylabel("Twist angle (^\circ)")
zlabel("Barrier (eV)")
set(gca, 'ydir', 'reverse')
view(45,60)