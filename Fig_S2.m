R = [100 150 200 250 300 400 500 600 700 800 1000];
angle1 = [7.46 6.12 5.31 4.72 4.32 3.79 3.38 3.09 2.83 2.68 2.38];
N_angle1 = [747 613 533 474 434 379 330 310 284 628 797];

f=figure;
f.Position = [0,0,600,600];
t=tiledlayout(4,4,"Padding","compact",'TileSpacing','tight');
nexttile([1,1]); hold on;
radius = 3e4:100:7e4;
angle = 74.2208./sqrt(radius)+0.0485;
plot(angle,radius,'r','LineWidth',1)
scatter(74.2208./sqrt(5e4)+0.0485,5e4,10,'ro','filled')
xlim([0,1.9])
xticks([])
yticks([3e4,5e4,7e4])
ylim([3e4,7e4])
ax = gca;
ax.XColor = "white";

nexttile([1,3]); hold on
x = linspace(0,1,1001);
y = cos(2*pi*9*x).*sin(x*pi);
plot(x,y,'k','LineWidth',1)
plot(x,sin(x*pi),':','color',[0.5,0.5,0.5])
plot(x,-sin(x*pi),':','color',[0.5,0.5,0.5])
xlim([-0.05,1.05])
ylim([-1.5,1.5])
xticks([])
yticks([])


nexttile([3,1]); hold on;
for i = 1:11
    r = R(i);
    E_MD = readmatrix("energy_R"+num2str(r)+".dat");
    E_MD(:,2) = (E_MD(:,2) - E_MD(N_angle1(i),2))/20 + r;
    plot(E_MD(:,1),E_MD(:,2),'-k','LineWidth',0.5)
end
xlim([0,1.9])
ylim([0,1200])


nexttile([3,3]); hold on;
for i = 1:11
    r = R(i);
    E_MD = readmatrix("energy_R"+num2str(r)+".dat");
    E_MD(:,2) = E_MD(:,2) - E_MD(N_angle1(i),2) + r;
    plot(E_MD(:,1),E_MD(:,2),'-k','LineWidth',0.5)
end

angle = 2:0.01:8;
radius = 74.2208^2./(angle-0.0485).^2;
plot(angle,radius,'r','LineWidth',1)

scatter(angle1,R,10,'ro','filled')

maxline = 1300*ones(size(angle));
inBetween = [radius, maxline];
x2 = [angle, fliplr(angle)];
fill(x2, inBetween, 'r','FaceAlpha',0.1,'EdgeColor','none');

xlim([2,8])
ylim([0,1200])
yticks([])
xlabel(t,"Twist angle ($^\circ$)",'interpreter','latex','fontsize',13.2)
ylabel(t,"Flake size $R/a_0$ $|$ Energy (eV)",'interpreter','latex','fontsize',13.2)


