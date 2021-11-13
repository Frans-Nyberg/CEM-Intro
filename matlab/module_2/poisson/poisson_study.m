savedir = 'results_2_1/';
save_image = @(filename) saveas(gca, strcat(savedir, filename), 'png');
args = struct();
nomial_order = @(hh,p,ref,sc) (hh/hh(end)).^p*ref*sc;

%% initial values, following Ch 3 in T.Rylander 2013
a = 0.01;   b = 0.01;   c = 0.02;   d = 0.02;
tol = 1e-9;
rel = 1.9;
%rel = 1;

%% 1.0 plot of potential
n = 20;
args.visualise = 1;
capacitor(a,b,c,d,n,tol,rel,args);
save_image("potential-GS")

%% 1.1, 1.2 convergence test
n = 10*15;
args.visualise = 0;
C15 = capacitor(a,b,c,d,n,tol,rel,args);

%%
N = 5;
nN = (10:10:10*N);
hN = 0.5*c./nN;
CN = zeros(N,1);
for k=1:N
    n = nN(k);
    CN(k) = capacitor(a,b,c,d,n,tol,rel,args);
end
%%
errN = abs(CN-C15);

%%
subplot(1,2,1)
plot(2*nN,CN,'-o')
xlabel('grid points in x-direction')
ylabel('capacitance')
subplot(1,2,2)
hold on
plot(log10(hN/2),log10(errN),'-o')
plot(log10(hN/2),log10(nomial_order(hN,1,errN(N),1)))
plot(log10(hN/2),log10(nomial_order(hN,2,errN(N),1)))
plot(log10(hN/2),log10(nomial_order(hN,3,errN(N),1)))
legend('C_{err}', 'O(h)', 'O(h^2)', 'O(h^3)')
hold off
xlabel('log10(h_x)')
ylabel('log10(error)')

%%
save_image("conv-GS")