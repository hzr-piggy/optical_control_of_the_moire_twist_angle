R = [300:100:800 1000];
const = [-70.204311497695898 -70.109038112391573 -70.050632045944724 ...
    -70.147107756162384 -70.100452568628896 -70.060299409742115 -70.016388919432387];
coeff = mean([0.395719587897704 0.395852440006244 0.394288238292406 ...
    0.395096462362051 0.395021012828176 0.393361450672660 0.394436252754282]);
x = linspace(0,10,101);
figure; hold on
for i = 1:7
    E_MD = readmatrix("energy_R"+num2str(R(i))+".dat");
    theta = E_MD(:,1);
    energy = (E_MD(:,2) - E_MD(1,2))/R(i)^2*1000;
    h(i) = scatter(theta,energy,'.','MarkerEdgeColor',(8-[i,i,i])/8,'MarkerFaceColor',(8-[i,i,i])/8);
    plot(x,const(i)+coeff*cosd(6*x),'color',(8-[i,i,i])/8,'LineWidth',2)
end
ylim([-0.0701,-0.0695]*1000)
legend(h,{'N = 300','N = 400','N = 500','N = 600','N = 700','N = 800','N = 1000'})
xlabel("Twist angle ($^\circ$)",'Interpreter','latex')
ylabel("Energy/N$^2$ (meV)",'Interpreter','latex')