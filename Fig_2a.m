clearvars;
f = figure('DefaultAxesFontSize',14,'Position',[0,0,500,400]);
% t=tiledlayout(4,4,"TileSpacing","compact","Padding","compact");
t=tiledlayout(3,4,"Padding","loose","TileSpacing","compact");

%% Variables 
chi = 1; % dielectric constant
E0 = 1; % peak electric field
ell = 1; % OAM number, choose 1 or -1, higher ell not supported
sigma = 0; % 0 for linear, +1/-1 for circular polarized light
ratio = 0.5; % ratio of beam radius to flake radius
norm_factor = E0*sqrt(4/pi);

%% Calculate electric field and polarization vectors on small flake
RADIUS = 5; % ratio of flake radius to lattice constant, unitless
w0 = RADIUS*ratio;

%% Build a hexagonal flake
coord_hex = get_flake(RADIUS);
a=RADIUS; b=a; c=a;
hexagon_x=[a c/2 -b/2 -a -c/2 b/2 a]; 
hexagon_y=[0 c*sqrt(3)/2 b*sqrt(3)/2 0 -c*sqrt(3)/2 -b*sqrt(3)/2 0];
lattice_vec = [1 0; -1/2 sqrt(3)/2];
coord_cart = coord_hex*lattice_vec; 
[coord_phi, coord_R]  = cart2pol(coord_cart(:,1),coord_cart(:,2));
[coord_x,coord_y] = pol2cart(coord_phi,coord_R);

%% Plot electric field and polarization vectors on a coarse grid
for i = 1:4
    time = (i-1)/4*2*pi;
    E = norm_factor*(coord_R/w0).*exp(-coord_R.^2/w0^2);
    E_x = E.*cos(-time-ell*coord_phi);
    E_y = sigma*E.*sin(-time-ell*coord_phi);
    index1 = E_x>=0;
    index2 = E_x<0;
    nexttile(i); hold on
    quiver(coord_x(index1),coord_y(index1),E_x(index1),E_y(index1),0.5,...
        'color',[1,0,0],'LineWidth',0.5,'Alignment','center');
    quiver(coord_x(index2),coord_y(index2),E_x(index2),E_y(index2),0.5,...
        'color',[0,0,1],'LineWidth',0.5,'Alignment','center');
    xticks([]); yticks([])
    xlim([-2.1*w0,2.1*w0]); ylim([-2.1*w0,2.1*w0]);
    plot(hexagon_x,hexagon_y,'color',[0.5,0.5,0.5]);
    axis off

    PIR_x = norm_factor*chi ...
        *(coord_R/w0).*exp(-coord_R.^2/w0^2)...
        .*cos(-(time+pi/2)-ell*coord_phi);
    PIR_y = sigma*norm_factor*chi ...
        *(coord_R/w0).*exp(-coord_R.^2/w0^2)...
        .*sin(-(time+pi/2)-ell*coord_phi);
    index1 = PIR_x>=0;
    index2 = PIR_x<0;
    nexttile(i+4); hold on
    quiver(coord_x(index1),coord_y(index1),PIR_x(index1),PIR_y(index1),0.5,...
        'color',[1,0,0],'LineWidth',0.5,'Alignment','center');
    quiver(coord_x(index2),coord_y(index2),PIR_x(index2),PIR_y(index2),0.5,...
        'color',[0,0,1],'LineWidth',0.5,'Alignment','center');
    xticks([]); yticks([])
    xlim([-2.1*w0,2.1*w0]); ylim([-2.1*w0,2.1*w0]);
    plot(hexagon_x,hexagon_y,'color',[0.5,0.5,0.5]);
    axis off
end
%% Calculate torque on a larger flake (finer grid)
RADIUS = 20; % ratio of flake radius to lattice constant, unitless
w0 = RADIUS*ratio;
coord_hex = get_flake(RADIUS);

a=RADIUS; b=a; c=a;
hexagon_x=[a c/2 -b/2 -a -c/2 b/2 a]; 
hexagon_y=[0 c*sqrt(3)/2 b*sqrt(3)/2 0 -c*sqrt(3)/2 -b*sqrt(3)/2 0];
lattice_vec = [1 0; -1/2 sqrt(3)/2];
coord_cart = coord_hex*lattice_vec; % unit = \AA
[coord_phi, coord_R]  = cart2pol(coord_cart(:,1),coord_cart(:,2));
[coord_x,coord_y] = pol2cart(coord_phi,coord_R);

for i = 1:4
    time = (i-1)/4*2*pi;
    E = norm_factor*(coord_R/w0).*exp(-coord_R.^2/w0^2);
    E_x = E.*cos(-time-ell*coord_phi);
    E_y = sigma*E.*sin(-time-ell*coord_phi);

    PIR_x = norm_factor*chi ...
        *(coord_R/w0).*exp(-coord_R.^2/w0^2)...
        .*cos(-(time+pi/2)-ell*coord_phi);
    PIR_y = sigma*norm_factor*chi ...
        *(coord_R/w0).*exp(-coord_R.^2/w0^2)...
        .*sin(-(time+pi/2)-ell*coord_phi);

    %% Charge contribution to torque
    rho_b = -norm_factor*chi*exp(-coord_R.^2/w0^2)...
        .*(cos(-(time+pi/2)-ell*coord_phi).*(1/w0-2*coord_R.^2/w0^3).*cos(coord_phi)... % x
        +cos(-(time+pi/2)+ell*pi/2-ell*coord_phi)*(1/w0).*sin(coord_phi)... % x
        +sigma*sin(-(time+pi/2)-ell*coord_phi).*(1/w0-2*coord_R.^2/w0^3).*sin(coord_phi)... % y
        -sigma*sin(-(time+pi/2)+ell*pi/2-ell*coord_phi)*(1/w0).*cos(coord_phi)); % y
    torque_E = (coord_x.*E_y-coord_y.*E_x).*rho_b;
    %% Current contribution to torque
    J_x = norm_factor*chi ...
        *(coord_R/w0).*exp(-coord_R.^2/w0^2)...
        .*cos(-(time+pi/2)-pi/2-ell*coord_phi);
    J_y = sigma*norm_factor*chi ...
        *(coord_R/w0).*exp(-coord_R.^2/w0^2)...
        .*sin(-(time+pi/2)-pi/2-ell*coord_phi);
    Bz = norm_factor*exp(-coord_R.^2/w0^2)...
        .*(cos(-time+pi/2-ell*coord_phi).*(1/w0-2*coord_R.^2/w0^3).*sin(coord_phi)...
        +cos(-time+pi/2-ell*pi/2-ell*coord_phi)*(1/w0).*cos(coord_phi) ...
        -sigma*sin(-time+pi/2-ell*coord_phi).*(1/w0-2*coord_R.^2/w0^3).*cos(coord_phi) ...
        +sigma*sin(-time+pi/2-ell*pi/2-ell*coord_phi)*(1/w0).*sin(coord_phi));
    torque_B = -Bz.*(coord_x.*J_x+coord_y.*J_y);

    %% Plot total torque
    torque_total = torque_E+torque_B;
    torque_max = max(abs(torque_total));
    colormap = [1 1 1]+(torque_total>0).*torque_total/torque_max*[0,-1,-1] ...
        +(torque_total<0).*torque_total/torque_max*[1,1,0];
    nexttile(i+8); hold on
    h=scatter3(coord_x,coord_y,torque_total,7,colormap,'filled');
    plot(hexagon_x,hexagon_y,'color',[0.5,0.5,0.5]);
    xticks([]); yticks([])
    xlim([-2.1*w0,2.1*w0]); ylim([-2.1*w0,2.1*w0]);
    axis off
end
disp("Maximum torque is "+num2str(torque_max))

% Function to construct a hexagonal flake
function coord = get_flake(RADIUS)
    N_CELL = 1+RADIUS*(RADIUS+1)/2*6;
    coord = zeros(N_CELL,2);
    start = 1;
    for R = 1:RADIUS
        for n = 1:R
            coord(start+n,:) = [R, n-1];
            coord((start+R+n),:) = [(R-n+1), R];
            coord((start+2*R+n),:) = [-(n-1),(R-n+1)];
            coord((start+3*R+n),:) = [-R, -(n-1)];
            coord((start+4*R+n),:) = [-(R-n+1), -R];
            coord((start+5*R+n),:) = [(n-1), -(R-n+1)];
        end
        start = start + 6*R;
    end
end
