f = figure('DefaultAxesFontSize',14,'Position',[0,50,600,270]);
t=tiledlayout(1,2,"TileSpacing","tight","Padding","compact");

p = 0; w0 = 0.5; l = 1;
nexttile; hold on;
W0 = 0.5;
r = linspace(0,1.1,101);
phi = linspace(0,2*pi,101);
[temp1,temp2] = meshgrid(r,phi);
[X,Y] = pol2cart(temp2,temp1);
I = (2*factorial(p)/pi/factorial(p+abs(l)))...
        *(2*(X.^2+Y.^2)/w0^2).^abs(l)...
        .*(laguerreL(p,abs(l),2*(X.^2+Y.^2)/w0^2)).^2 ...
        .*exp(-2*(X.^2+Y.^2)/w0^2);
contourf(X,Y,I,30,'EdgeColor','none');
a = 1; b = 1; c = 1;
hexagon_x=[a c/2 -b/2 -a -c/2 b/2 a]; hexagon_y=[0 c*sqrt(3)/2 b*sqrt(3)/2 0 -c*sqrt(3)/2 -b*sqrt(3)/2 0];
plot(hexagon_x,hexagon_y,'color',[0.5,0.5,0.5])
xlim([-2.2*W0,2.2*W0]); ylim([-2.2*W0,2.2*W0]);
axis off
c=colorbar('eastoutside','xticklabel',["0","max"],'xtick',[0,max(max(I))-0.01]);

nexttile; hold on
x = linspace(0,2,101);
style = ["-","--",":"];
for l = 1:3
    Lx = (2*factorial(p)/pi/factorial(p+abs(l)))...
        *(sqrt(2)*x/w0).^abs(2*l)...
        .*(laguerreL(p,abs(l),2*x.^2/w0^2)).^2 ...
        .*exp(-2*x.^2/w0^2);
    plot(x,Lx,'k','LineWidth',2,'LineStyle',style(l))
end

xticks([0,0.5,1])
xlim([0,1.25])
xticklabels({"0","w_0","R"})
yticks([])
ylabel("Torque")
legend({"$\ell$=1","$\ell$=2","$\ell$=3"},'interpreter','latex')