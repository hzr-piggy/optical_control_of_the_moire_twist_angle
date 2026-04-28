R = [100 150 200 250 300 400 500 600 700 800 1000];
angle1 = [7.46 6.12 5.31 4.72 4.32 3.79 3.38 3.09 2.83 2.68 2.38];
N_angle1 = [747 613 533 474 434 379 330 310 284 628 797];

f=figure;
f.Position = [0,0,800,800];
hold on;
for i = 1:11
    r = R(i);
    E_MD = readmatrix("energy_R"+num2str(r)+".dat");
    E_MD(:,2) = E_MD(:,2) - E_MD(N_angle1(i),2) + r;
    plot(E_MD(:,1),E_MD(:,2),'-k','LineWidth',0.5)
end

angle = 2:0.01:10;
radius = 74.2208^2./(angle-0.0485).^2;
plot(angle,radius,'r','LineWidth',1)

scatter(angle1,R,10,'ro','filled')

coeff = mean([0.395719587897704 0.395852440006244 0.394288238292406 ...
    0.395096462362051 0.395021012828176 0.393361450672660 0.394436252754282])/1000;
shift = [0.8390 1.6939 2.4146 3.0522 3.7646 4.9766 8.6175 7.5285 8.7879 10.7424 11.4477];
theta = linspace(2,10,10001);
for i = 1:11
    r = R(i);
    E_fit = 1439/2*theta.^(-3).*cosd(0.1103*r*theta*180/pi)+coeff*r^2*(cosd(theta*6)-1)+r+shift(i);
    plot(theta,E_fit,'-','LineWidth',0.5,'Color', '#FFA500')
end

h(1) = plot([0,1],[-1,-1],'-k','LineWidth',0.5);
h(2) = plot([0,1],[-1,-1],'-c','LineWidth',0.5,'Color', '#FFA500');
legend(h,["From simulation","From Eq. S4"])

xlim([2,10])
ylim([0,1100])
xlabel("Twist angle ($^\circ$)",'interpreter','latex','fontsize',15)
ylabel("Flake size $R/a_0$ $|$ Energy (eV)",'interpreter','latex','fontsize',15)


