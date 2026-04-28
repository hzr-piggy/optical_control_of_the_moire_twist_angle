dt = [0.1,0.5,1]*2*pi;
f = figure; 
tile = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
opts = odeset('RelTol',1e-6);
nexttile; hold on
Ji = 0.75; 
rep_max = [2,6,2];
extend = [true,false,false];
markers = ['o','x','+'];
style = ["-","--",":"];
for i = 1:3
    dti = dt(i); 
    rep = 0;
    initial = [0,2*Ji];
    w = [0,0];
    while rep<rep_max(i)
        [t,w] = ode45(@(t,w) derivatives(t,w,0.0),...
                    [0,dti],initial,opts);
        initial = w(end,:) + [0, 2*Ji];
        plot((t+rep*dti)/2/pi,w(:,1),'k','LineWidth',1,'LineStyle',style(i));
        scatter(rep*dti/2/pi,w(1,1),20,'red','Marker',markers(i),'LineWidth',1)
        rep = rep + 1;
    end
    if extend(i)
        initial = w(end,:) ;
        [t,w] = ode45(@(t,w) derivatives(t,w,0.0),...
                            [0,0.4*pi],initial,opts);
        plot((t+rep*dti)/2/pi,w(:,1),'k','LineWidth',1,'LineStyle',style(i));
    end
end
ylim([-2,4])
yticks([-pi/2,0,pi/2,pi])
yticklabels({"-\pi/2","0","\pi/2","\pi"})
yline(0,'-')
yline(pi,'-')
xticks([])
title("(i)",'fontsize',14)


nexttile; hold on
J = [0.75,0.35,0.1];
rep_max = [2,4,30];
extend = [true,true,false];
markers = ['o','s','*'];
dti = 0.1*2*pi;
for i = 1:3
    Ji = J(i); 
    rep = 0;
    initial = [0,2*Ji];
    w = [0,0];
    while rep<rep_max(i)
        [t,w] = ode45(@(t,w) derivatives(t,w,0.0),...
                    [0,dti],initial,opts);
        initial = w(end,:) + [0, 2*Ji];
        plot((t+rep*dti)/2/pi,w(:,1),'k','LineWidth',1,'LineStyle',style(i));
        scatter(rep*dti/2/pi,w(1,1),20,'red','Marker',markers(i),'LineWidth',1)
        rep = rep + 1;
    end
    if extend(i)
        initial = w(end,:) ;
        [t,w] = ode45(@(t,w) derivatives(t,w,0.0),...
                            [0,0.6*pi],initial,opts);
        plot((t+rep*dti)/2/pi,w(:,1),'k','LineWidth',1,'LineStyle',style(i));
    end
end
title("(ii)",'fontsize',14)
ylim([-2,4])
yticks([-pi/2,0,pi/2,pi])
yticklabels({"-\pi/2","0","\pi/2","\pi"})
yline(0,'-')
yline(pi,'-')
xlabel(tile,"Time $t\ (T_0)$",'interpreter','latex','fontsize',14)
ylabel(tile,"Rescaled angle $\Theta'=k\Theta$",'interpreter','latex','fontsize',14)

function derivs = derivatives(tf,wf,GAMMA)
    theta = wf(1); % theta = rotation angle in degree
    theta_dot = wf(2);
    theta_2dot = - sin(theta) - GAMMA*theta_dot;
    derivs = [theta_dot; theta_2dot];
end