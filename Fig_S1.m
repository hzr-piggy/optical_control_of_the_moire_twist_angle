clearvars;
f = figure('DefaultAxesFontSize',14,'Position',[0,0,400,800]);
t=tiledlayout(8,4,"Padding","compact","TileSpacing","tight");

sigma = -1; % spin angular momentum
% sigma takes value of 0 (linear) and +/-1 (circular polarized)
chi = 1; % dielectric susceptibility
E0 = 1; % strengh of peak electric field
omegaIR = 1; % frequency of IR active phonon
ell = 1; % orbital angular momentum index
torque_max = 0.2342; % used to scale the colorbar
if sigma == 0
    norm_factor = E0*sqrt(4/pi);
else
    norm_factor = E0*sqrt(4/pi)/sqrt(2);
end

%% Plot electric field and polarization, size = 10
RADIUS = 5; % ratio of flake radius to lattice constant, unitless
ratio = 0.5;
LATTICE_CONST = 2.5; % \AA
N_CELL = 1+RADIUS*(RADIUS+1)/2*6;
AREA_CELL = LATTICE_CONST^2/2*sqrt(3); % \AA^2 
W0 = LATTICE_CONST*RADIUS*ratio;
coord_hex = zeros(N_CELL,2);
start = 1;
for R = 1:RADIUS
    for n = 1:R
        coord_hex(start+n,:) = [R, n-1];
        coord_hex((start+R+n),:) = [(R-n+1), R];
        coord_hex((start+2*R+n),:) = [-(n-1),(R-n+1)];
        coord_hex((start+3*R+n),:) = [-R, -(n-1)];
        coord_hex((start+4*R+n),:) = [-(R-n+1), -R];
        coord_hex((start+5*R+n),:) = [(n-1), -(R-n+1)];
    end
    start = start + 6*R;
end
a=LATTICE_CONST*RADIUS; b=a; c=a;
hexagon_x=[a c/2 -b/2 -a -c/2 b/2 a]; hexagon_y=[0 c*sqrt(3)/2 b*sqrt(3)/2 0 -c*sqrt(3)/2 -b*sqrt(3)/2 0];
lattice_vec = [1 0; -1/2 sqrt(3)/2]*LATTICE_CONST;
COORD_CART = coord_hex*lattice_vec; % unit = \AA
[COORD_PHI, COORD_R]  = cart2pol(COORD_CART(:,1),COORD_CART(:,2));

% Get surface and edge coordinates and electric field vectors
[coord_x,coord_y] = pol2cart(COORD_PHI,COORD_R);
n_edge = RADIUS*6; n_surf = N_CELL - n_edge;
coord_surf_x = COORD_CART(1:n_surf,1);
coord_surf_y = COORD_CART(1:n_surf,2);
coord_edge_x = COORD_CART((n_surf+1):end,1);
coord_edge_y = COORD_CART((n_surf+1):end,2);
edges = [sqrt(3)/2 1/2; 0 1; -sqrt(3)/2 1/2; ...
    -sqrt(3)/2 -1/2; 0 -1; sqrt(3)/2 -1/2];
corners = [1 0; 1/2 sqrt(3)/2; -1/2 sqrt(3)/2; ...
    -1 0; -1/2 -sqrt(3)/2; 1/2 -sqrt(3)/2];
edge_vec_x = zeros(6*RADIUS,1);
edge_vec_y = zeros(6*RADIUS,1);
edge_vec_x(1:RADIUS:(6*RADIUS)) = corners(:,1);
edge_vec_y(1:RADIUS:(6*RADIUS)) = corners(:,2);
for edge = 0:5
    edge_vec_x((edge*RADIUS+2):(edge*RADIUS+RADIUS)) = edges(edge+1,1);
    edge_vec_y((edge*RADIUS+2):(edge*RADIUS+RADIUS)) = edges(edge+1,2);
end

for i = 1:4
    time = (i-1)/4*2*pi;
    E = norm_factor*(COORD_R/W0).*exp(-COORD_R.^2/W0^2);
    E_x = E.*cos(-time-ell*COORD_PHI);
    E_y = sigma*E.*sin(-time-ell*COORD_PHI);
    nexttile(i); hold on
    quiver(coord_x,coord_y,E_x,E_y,0.5,'color',[0,0,0],'LineWidth',0.5,'Alignment','center');
    xticks([]); yticks([])
    xlim([-2.1*W0,2.1*W0]); ylim([-2.1*W0,2.1*W0]);
    plot(hexagon_x,hexagon_y,'color',[0.5,0.5,0.5]);
    axis off

    PIR_x = norm_factor*chi ...
        *(COORD_R/W0).*exp(-COORD_R.^2/W0^2)...
        .*cos(-(time+pi/2)-ell*COORD_PHI);
    PIR_y = sigma*norm_factor*chi ...
        *(COORD_R/W0).*exp(-COORD_R.^2/W0^2)...
        .*sin(-(time+pi/2)-ell*COORD_PHI);
    
    J_x = norm_factor*chi ...
        *(COORD_R/W0).*exp(-COORD_R.^2/W0^2)...
        .*cos(-(time+pi/2)-pi/2-ell*COORD_PHI);
    J_y = sigma*norm_factor*chi ...
        *(COORD_R/W0).*exp(-COORD_R.^2/W0^2)...
        .*sin(-(time+pi/2)-pi/2-ell*COORD_PHI);
    J_surf_x = J_x(1:n_surf);
    J_surf_y = J_y(1:n_surf);
    nexttile(i+12); hold on
    quiver(coord_x,coord_y,J_x,J_y,0.5,'color',[0,0,0],'LineWidth',0.5,'Alignment','center');
    xticks([]); yticks([])
    xlim([-2.1*W0,2.1*W0]); ylim([-2.1*W0,2.1*W0]);
    plot(hexagon_x,hexagon_y,'color',[0.5,0.5,0.5]);
    axis off

end
%% Plot torque, size = 20
RADIUS = 20; % ratio of flake radius to lattice constant, unitless
N_CELL = 1+RADIUS*(RADIUS+1)/2*6;
AREA_CELL = LATTICE_CONST^2/2*sqrt(3); % \AA^2 
W0 = LATTICE_CONST*RADIUS*ratio;
coord_hex = zeros(N_CELL,2);
start = 1;
for R = 1:RADIUS
    for n = 1:R
        coord_hex(start+n,:) = [R, n-1];
        coord_hex((start+R+n),:) = [(R-n+1), R];
        coord_hex((start+2*R+n),:) = [-(n-1),(R-n+1)];
        coord_hex((start+3*R+n),:) = [-R, -(n-1)];
        coord_hex((start+4*R+n),:) = [-(R-n+1), -R];
        coord_hex((start+5*R+n),:) = [(n-1), -(R-n+1)];
    end
    start = start + 6*R;
end
a=LATTICE_CONST*RADIUS; b=a; c=a;
hexagon_x=[a c/2 -b/2 -a -c/2 b/2 a]; hexagon_y=[0 c*sqrt(3)/2 b*sqrt(3)/2 0 -c*sqrt(3)/2 -b*sqrt(3)/2 0];
lattice_vec = [1 0; -1/2 sqrt(3)/2]*LATTICE_CONST;
COORD_CART = coord_hex*lattice_vec; % unit = \AA
[COORD_PHI, COORD_R]  = cart2pol(COORD_CART(:,1),COORD_CART(:,2));

% Get surface and edge coordinates and electric field vectors
[coord_x,coord_y] = pol2cart(COORD_PHI,COORD_R);
n_edge = RADIUS*6; n_surf = N_CELL - n_edge;
coord_surf_x = COORD_CART(1:n_surf,1);
coord_surf_y = COORD_CART(1:n_surf,2);
coord_edge_x = COORD_CART((n_surf+1):end,1);
coord_edge_y = COORD_CART((n_surf+1):end,2);
edges = [sqrt(3)/2 1/2; 0 1; -sqrt(3)/2 1/2; ...
    -sqrt(3)/2 -1/2; 0 -1; sqrt(3)/2 -1/2];
corners = [1 0; 1/2 sqrt(3)/2; -1/2 sqrt(3)/2; ...
    -1 0; -1/2 -sqrt(3)/2; 1/2 -sqrt(3)/2];
edge_vec_x = zeros(6*RADIUS,1);
edge_vec_y = zeros(6*RADIUS,1);
edge_vec_x(1:RADIUS:(6*RADIUS)) = corners(:,1);
edge_vec_y(1:RADIUS:(6*RADIUS)) = corners(:,2);
for edge = 0:5
    edge_vec_x((edge*RADIUS+2):(edge*RADIUS+RADIUS)) = edges(edge+1,1);
    edge_vec_y((edge*RADIUS+2):(edge*RADIUS+RADIUS)) = edges(edge+1,2);
end

for i = 1:4
    time = (i-1)/4*2*pi;
    E = norm_factor*(COORD_R/W0).*exp(-COORD_R.^2/W0^2);
    E_x = E.*cos(-time-ell*COORD_PHI);
    E_y = sigma*E.*sin(-time-ell*COORD_PHI);
    E_surf_x = E_x(1:n_surf);
    E_surf_y = E_y(1:n_surf);
    E_edge_x = E_x((n_surf+1):end);
    E_edge_y = E_y((n_surf+1):end);

    PIR_x = norm_factor*chi ...
        *(COORD_R/W0).*exp(-COORD_R.^2/W0^2)...
        .*cos(-(time+pi/2)-ell*COORD_PHI);
    PIR_y = sigma*norm_factor*chi ...
        *(COORD_R/W0).*exp(-COORD_R.^2/W0^2)...
        .*sin(-(time+pi/2)-ell*COORD_PHI);
    PIR_surf_x = PIR_x(1:n_surf);
    PIR_surf_y = PIR_y(1:n_surf);
    PIR_edge_x = PIR_x((n_surf+1):end);
    PIR_edge_y = PIR_y((n_surf+1):end);

    %% Charge contribution
    rho_b = -norm_factor*chi*exp(-COORD_R.^2/W0^2)...
        .*(cos(-(time+pi/2)-ell*COORD_PHI).*(1/W0-2*COORD_R.^2/W0^3).*cos(COORD_PHI)... % x
        +cos(-(time+pi/2)+ell*pi/2-ell*COORD_PHI)*(1/W0).*sin(COORD_PHI)... % x
        +sigma*sin(-(time+pi/2)-ell*COORD_PHI).*(1/W0-2*COORD_R.^2/W0^3).*sin(COORD_PHI)... % y
        -sigma*sin(-(time+pi/2)+ell*pi/2-ell*COORD_PHI)*(1/W0).*cos(COORD_PHI)); % y
    rho_b = rho_b(1:n_surf);
    nexttile(i+4); hold on
    rho_max = 0.0638;
    colormap_rho = [1 1 1]+(rho_b>0).*rho_b/rho_max*[0,-1,-1] ...
        +(rho_b<0).*rho_b/rho_max*[1,1,0];
    scatter3(coord_surf_x,coord_surf_y,rho_b,7,colormap_rho,'filled')
    xticks([]); yticks([])
    xlim([-2.1*W0,2.1*W0]); ylim([-2.1*W0,2.1*W0]);
    axis off

    % Compute bound charges on the edge
    sigma_b = PIR_edge_x.*edge_vec_x...
        +PIR_edge_y.*edge_vec_y; % unit e/AA

    torque_E = (coord_surf_x.*E_surf_y-coord_surf_y.*E_surf_x).*rho_b;
    colormapE = [1 1 1]+(torque_E>0).*torque_E/torque_max*[0,-1,-1] ...
        +(torque_E<0).*torque_E/torque_max*[1,1,0];
    torque_E_edge = (coord_edge_x.*E_edge_y-coord_edge_y.*E_edge_x).*sigma_b*LATTICE_CONST/AREA_CELL;
    colormapE_edge = [1 1 1]+(torque_E_edge>0).*torque_E_edge/torque_max*[0,-1,-1] ...
        +(torque_E_edge<0).*torque_E_edge/torque_max*[1,1,0];
    nexttile(i+8); hold on
    scatter3(coord_surf_x,coord_surf_y,torque_E,7,colormapE,'filled')
    scatter3(coord_edge_x,coord_edge_y,torque_E_edge,7,colormapE_edge,'filled')
    xticks([]); yticks([])
    xlim([-2.1*W0,2.1*W0]); ylim([-2.1*W0,2.1*W0]);
    axis off

    %% Current contribution
    J_x = norm_factor*chi ...
        *(COORD_R/W0).*exp(-COORD_R.^2/W0^2)...
        .*cos(-(time+pi/2)-pi/2-ell*COORD_PHI);
    J_y = sigma*norm_factor*chi ...
        *(COORD_R/W0).*exp(-COORD_R.^2/W0^2)...
        .*sin(-(time+pi/2)-pi/2-ell*COORD_PHI);
    J_surf_x = J_x(1:n_surf);
    J_surf_y = J_y(1:n_surf);

    Bz = norm_factor*exp(-COORD_R.^2/W0^2)...
        .*(cos(-time+pi/2-ell*COORD_PHI).*(1/W0-2*COORD_R.^2/W0^3).*sin(COORD_PHI)...
        +cos(-time+pi/2-ell*pi/2-ell*COORD_PHI)*(1/W0).*cos(COORD_PHI) ...
        -sigma*sin(-time+pi/2-ell*COORD_PHI).*(1/W0-2*COORD_R.^2/W0^3).*cos(COORD_PHI) ...
        +sigma*sin(-time+pi/2-ell*pi/2-ell*COORD_PHI)*(1/W0).*sin(COORD_PHI));
    Bz_max = 0.0638;
    Bz = Bz(1:n_surf);
    nexttile(i+16); hold on
    colormapBz = [1 1 1]+(Bz>0).*Bz/Bz_max*[0,-1,-1] ...
        +(Bz<0).*Bz/Bz_max*[1,1,0];
    scatter3(coord_surf_x,coord_surf_y,Bz,7,colormapBz,'filled')
    xticks([]); yticks([])
    xlim([-2.1*W0,2.1*W0]); ylim([-2.1*W0,2.1*W0]);
    axis off

    torque_B = -Bz.*(coord_surf_x.*J_surf_x+coord_surf_y.*J_surf_y);
    colormapB = [1 1 1]+(torque_B>0).*torque_B/torque_max*[0,-1,-1] ...
        +(torque_B<0).*torque_B/torque_max*[1,1,0];
    
    nexttile(i+20); hold on
    scatter3(coord_surf_x,coord_surf_y,torque_B,7,colormapB,'filled')
    xticks([]); yticks([])
    xlim([-2.1*W0,2.1*W0]); ylim([-2.1*W0,2.1*W0]);
    axis off

    %% Total torque
    torque_total = torque_E+torque_B;
    colormap = [1 1 1]+(torque_total>0).*torque_total/torque_max*[0,-1,-1] ...
        +(torque_total<0).*torque_total/torque_max*[1,1,0];
    nexttile(i+24); hold on
    scatter3(coord_surf_x,coord_surf_y,torque_total,7,colormap,'filled')
    xticks([]); yticks([])
    xlim([-2.1*W0,2.1*W0]); ylim([-2.1*W0,2.1*W0]);
    scatter3(coord_edge_x,coord_edge_y,torque_E_edge,7,colormapE_edge,'filled')
    axis off
end

RADIUS = 100; % use larger flake size for integration
N_CELL = 1+RADIUS*(RADIUS+1)/2*6;
AREA_CELL = LATTICE_CONST^2/2*sqrt(3); % \AA^2 
W0 = LATTICE_CONST*RADIUS*ratio;
coord_hex = zeros(N_CELL,2);
start = 1;
for R = 1:RADIUS
    for n = 1:R
        coord_hex(start+n,:) = [R, n-1];
        coord_hex((start+R+n),:) = [(R-n+1), R];
        coord_hex((start+2*R+n),:) = [-(n-1),(R-n+1)];
        coord_hex((start+3*R+n),:) = [-R, -(n-1)];
        coord_hex((start+4*R+n),:) = [-(R-n+1), -R];
        coord_hex((start+5*R+n),:) = [(n-1), -(R-n+1)];
    end
    start = start + 6*R;
end
a=LATTICE_CONST*RADIUS; b=a; c=a;
hexagon_x=[a c/2 -b/2 -a -c/2 b/2 a]; hexagon_y=[0 c*sqrt(3)/2 b*sqrt(3)/2 0 -c*sqrt(3)/2 -b*sqrt(3)/2 0];
lattice_vec = [1 0; -1/2 sqrt(3)/2]*LATTICE_CONST;
COORD_CART = coord_hex*lattice_vec; % unit = \AA
[COORD_PHI, COORD_R]  = cart2pol(COORD_CART(:,1),COORD_CART(:,2));

% Get surface and edge coordinates and electric field vectors
[coord_x,coord_y] = pol2cart(COORD_PHI,COORD_R);
n_edge = RADIUS*6; n_surf = N_CELL - n_edge;
coord_surf_x = COORD_CART(1:n_surf,1);
coord_surf_y = COORD_CART(1:n_surf,2);
coord_edge_x = COORD_CART((n_surf+1):end,1);
coord_edge_y = COORD_CART((n_surf+1):end,2);
edges = [sqrt(3)/2 1/2; 0 1; -sqrt(3)/2 1/2; ...
    -sqrt(3)/2 -1/2; 0 -1; sqrt(3)/2 -1/2];
corners = [1 0; 1/2 sqrt(3)/2; -1/2 sqrt(3)/2; ...
    -1 0; -1/2 -sqrt(3)/2; 1/2 -sqrt(3)/2];
edge_vec_x = zeros(6*RADIUS,1);
edge_vec_y = zeros(6*RADIUS,1);
edge_vec_x(1:RADIUS:(6*RADIUS)) = corners(:,1);
edge_vec_y(1:RADIUS:(6*RADIUS)) = corners(:,2);
for edge = 0:5
    edge_vec_x((edge*RADIUS+2):(edge*RADIUS+RADIUS)) = edges(edge+1,1);
    edge_vec_y((edge*RADIUS+2):(edge*RADIUS+RADIUS)) = edges(edge+1,2);
end


Time = linspace(-1/8,7/8,101);
Torque = zeros(101,3);
for i = 1:101
    time = Time(i)*2*pi;
    E = norm_factor*(COORD_R/W0).*exp(-COORD_R.^2/W0^2);
    E_x = E.*cos(-time-ell*COORD_PHI);
    E_y = sigma*E.*sin(-time-ell*COORD_PHI);
    E_surf_x = E_x(1:n_surf);
    E_surf_y = E_y(1:n_surf);
    E_edge_x = E_x((n_surf+1):end);
    E_edge_y = E_y((n_surf+1):end);

    PIR_x = norm_factor*chi ...
        *(COORD_R/W0).*exp(-COORD_R.^2/W0^2)...
        .*cos(-(time+pi/2)-ell*COORD_PHI);
    PIR_y = sigma*norm_factor*chi ...
        *(COORD_R/W0).*exp(-COORD_R.^2/W0^2)...
        .*sin(-(time+pi/2)-ell*COORD_PHI);
    PIR_surf_x = PIR_x(1:n_surf);
    PIR_surf_y = PIR_y(1:n_surf);
    PIR_edge_x = PIR_x((n_surf+1):end);
    PIR_edge_y = PIR_y((n_surf+1):end);

    %% Charge contribution
    rho_b = -norm_factor*chi*exp(-COORD_R.^2/W0^2)...
        .*(cos(-(time+pi/2)-ell*COORD_PHI).*(1/W0-2*COORD_R.^2/W0^3).*cos(COORD_PHI)... % x
        +cos(-(time+pi/2)+ell*pi/2-ell*COORD_PHI)*(1/W0).*sin(COORD_PHI)... % x
        +sigma*sin(-(time+pi/2)-ell*COORD_PHI).*(1/W0-2*COORD_R.^2/W0^3).*sin(COORD_PHI)... % y
        -sigma*sin(-(time+pi/2)+ell*pi/2-ell*COORD_PHI)*(1/W0).*cos(COORD_PHI)); % y
    rho_b = rho_b(1:n_surf);

    % Compute bound charges on the edge
    sigma_b = PIR_edge_x.*edge_vec_x...
        +PIR_edge_y.*edge_vec_y; % unit e/AA

    torque_E = (coord_surf_x.*E_surf_y-coord_surf_y.*E_surf_x).*rho_b;
    torque_sigma = (coord_edge_x.*E_edge_y-coord_edge_y.*E_edge_x).*sigma_b/AREA_CELL*LATTICE_CONST;
    Torque(i,1) = sum(torque_E);
    Torque(i,2) = sum(torque_sigma);

    %% Current contribution
    J_x = norm_factor*chi ...
        *(COORD_R/W0).*exp(-COORD_R.^2/W0^2)...
        .*cos(-(time+pi/2)-pi/2-ell*COORD_PHI);
    J_y = sigma*norm_factor*chi ...
        *(COORD_R/W0).*exp(-COORD_R.^2/W0^2)...
        .*sin(-(time+pi/2)-pi/2-ell*COORD_PHI);
    J_surf_x = J_x(1:n_surf);
    J_surf_y = J_y(1:n_surf);

    Bz = norm_factor*exp(-COORD_R.^2/W0^2)...
        .*(cos(-time+pi/2-ell*COORD_PHI).*(1/W0-2*COORD_R.^2/W0^3).*sin(COORD_PHI)...
        +cos(-time+pi/2-ell*pi/2-ell*COORD_PHI)*(1/W0).*cos(COORD_PHI) ...
        -sigma*sin(-time+pi/2-ell*COORD_PHI).*(1/W0-2*COORD_R.^2/W0^3).*cos(COORD_PHI) ...
        +sigma*sin(-time+pi/2-ell*pi/2-ell*COORD_PHI)*(1/W0).*sin(COORD_PHI));

    torque_B = -Bz(1:n_surf).*(coord_surf_x.*J_surf_x+coord_surf_y.*J_surf_y);
    Torque(i,3) = sum(torque_B);
   
end
nexttile([1,4])
h(1)=plot(Time,Torque(:,1),'k--','LineWidth',1);
hold on
h(2)=plot(Time,Torque(:,3),'k:','LineWidth',1);
h(3)=plot(Time,Torque(:,2),'--','color',[0.5,0.5,0.5],'LineWidth',1);
h(4)=plot(Time,Torque(:,1)+Torque(:,2)+Torque(:,3),'k-','LineWidth',1);
yline(0)
ylim([-500,3500])
xlim([-1/8,7/8])
xticks([0,1/4,1/2,3/4])
xticklabels({"0","1/4 T","1/2 T","3/4 T"})
yticks(0)
legend(h,["Surf, charge","Surf, current", "Edge, charge","Total"])
%}