EAA = load("energy_R100.dat");
EAB = load("energy_R100_AB.dat");
figure('Position',[100,100,800,400]);  
tiledlayout(1,2,"TileSpacing","tight","Padding","tight")
nexttile; hold on
plot(EAA(:,1),EAA(:,2)-EAA(1,2),'k','LineWidth',1)
plot(EAB(:,1),EAB(:,2)-EAA(1,2),'r','LineWidth',1)
xlim([0,5])
ylabel("Energy (eV)",'Interpreter','latex')
xlabel("Twist angle ($^\circ$)",'Interpreter','latex')
legend(["AA-stacking","AB-stacking"])
nexttile; hold on
plot(EAA(:,1),EAA(:,2)-EAA(1,2),'k','LineWidth',1)
plot(EAB(:,1),EAB(:,2)-EAA(1,2),'r','LineWidth',1)
xlim([5,10])
xlabel("Twist angle ($^\circ$)",'Interpreter','latex')