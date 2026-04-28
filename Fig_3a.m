%% Without damping
% Plot number of pulses required to achieve twist
% as a function of interpulse period and effective angular momentum change.
%{
map = repmat((linspace(0.2,1,10))',1,3);
na = 240; nb = 500; maxrep = 10;
j = linspace(0.01,2.4,na)/2; % rescale by 2 according to Eq. (6)
dtau = linspace(0.02,10,nb)/2/pi;
[J,DTAU] = meshgrid(j,dtau);
data1 = load("num_pulses_damp0.mat").num_pulses;
figure('Position',[100,100,600,400]);
contourf(DTAU,J,data1,'LineColor','none');
colormap("pink");
hold on
yline(0.75,'--')
xline(0.1,'--')
scatter(0.1,0.75,'ok','LineWidth',1)
scatter(0.1,0.35,'sk','LineWidth',1)
scatter(0.1,0.1,'*k','LineWidth',1)
scatter(0.5,0.75,'xk','LineWidth',1)
scatter(1,0.75,'+k','LineWidth',1)
ylim([0,1.1])

cbh = colorbar ; %Create Colorbar
cbh.Ticks = linspace(1.5,9.5,10) ; 
cbh.TickLabels = {"1","2","3","4","5","6","7","8","9",">10"};
cbh.Label.String = 'Min number of pulses required';
ylabel("Effective angular momentum $J'$",'interpreter','latex','fontsize',14)
xlabel("Interpulse period $\Delta t (T_0)$",'interpreter','latex','fontsize',14)


%% With damping, Supplement
%{
na = 240; nb = 500; maxrep = 10;
j = linspace(0.01,2.4,na)/2;
dtau = linspace(0.02,10,nb)/2/pi;
[J,DTAU] = meshgrid(j,dtau);
data1 = load("num_pulses_damp0.mat").num_pulses;
data2 = load("num_pulses_damp0.1.mat").num_pulses;
figure("Position",[0,0,800,400]);
t=tiledlayout(1,2,"TileSpacing","none");
nexttile
contourf(DTAU,J,data2,'LineColor','none')
colormap("pink")
set(gca, 'Xdir','reverse')

title("$\gamma=0.1\Omega_0$",'Interpreter','latex','FontSize',14)
nexttile
contourf(DTAU,J,data1,'LineColor','none');
hold on
plot([0,0.5],[0,1],'r:')
yticklabels({})
cbh = colorbar ; %Create Colorbar
cbh.Ticks = linspace(1.5,9.5,10) ; 
cbh.TickLabels = {"1","2","3","4","5","6","7","8","9",">10"};
cbh.Label.String = 'Min number of pulses required';
title("$\gamma=0$",'Interpreter','latex','FontSize',14)
ylabel(t,"Effective ngular momentum $J'$",'interpreter','latex','fontsize',14)
xlabel(t,"Interpulse period $\Delta t\ (T_0)$",'interpreter','latex','fontsize',14)

insetAxes = axes('Position', [0.57 0.2 0.15 0.15]); % [x y width height]
box on; % Add a box around the inset
x = linspace(0,pi,101);
plot(x,1-cos(x), '-k');
hold on
plot(x,2/pi*x, '--k');
xticks([]); yticks([]);
xlabel("x"); ylabel("Energy")
%}

%% Generate num_pulses_damp0.mat and num_pulses_damp0.1.mat
%{
na = 240; nb = 500; maxrep = 10;
a = linspace(0.01,2.4,na);
dtau = linspace(0.02,10,nb);
[A,DTAU] = meshgrid(a,dtau);
num_pulses = ones(nb,na)*maxrep;
damping = 0;
% Find min number of pulses required to overcome the barrier
opts = odeset('RelTol',1e-6);
parfor i = 1:nb
    for j = 1:na
        initial = [0,A(i,j)];
        rep = 1; 
        while rep <= maxrep
            [~,w_theta_long] = ode45(@(t,w) derivatives(t,w,damping),...
                    [0,10],initial,opts);
            if max(w_theta_long(:,1)) > pi
                num_pulses(i,j) = rep;
                break
            else
                [~,w_theta] = ode45(@(t,w) derivatives(t,w,damping),...
                    [0,DTAU(i,j)],initial,opts);
                initial = w_theta(end,:) + [0, A(i,j)];
                rep = rep + 1;
            end
        end
    end
end
save("num_pulses_damp0.mat", "num_pulses")

function derivs = derivatives(tf,wf,GAMMA)
    theta = wf(1); % theta = rotation angle in degree
    theta_dot = wf(2);
    theta_2dot = -GAMMA*theta_dot ... 
        - sin(theta);
    derivs = [theta_dot; theta_2dot];
end
%}