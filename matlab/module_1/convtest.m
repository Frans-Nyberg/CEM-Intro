savedir = 'results_1/';
save_image = @(filename) saveas(gca, strcat(savedir, filename), 'png');

%% Following ch 2 in T.Rylander 2013
z = 1;
a = 1;
rule = 'simpson';
Iexact = 0.79335912;

%% Estimate of exact solution
I13 = integr(z, a, 2^13, rule);

%% Convergence test
expn0 = 1;
N = 12 - expn0 + 1;
expn = 1:1:N;
nN = 2.^expn;
hN = a./nN;

I = zeros(length(N),1);
for k=1:N
    n = nN(k);
    I(k) = integr(z, a, n, rule);
end

errN = abs(I-I13);

%% Extrapolation
Ifit = I(1:end); 
px = hN(1:end).^4;    % simpson
% try various orders
mM = 1:1:4;
M = length(mM);
assert(N>M)
IextM = zeros(M,1);
for k=1:M
    m = mM(k);
    pfitm = polyfit(px, Ifit, m);
    IextM(k) = pfitm(end);
end

%% Plots
grid on

figure(2)
subplot(2,1,1)
plot(mM, IextM, '-o')
xlabel('order')
ylabel('matlab cuts this out!')
title('zero of extrapolation')
subplot(2,1,2)
semilogy(mM, abs(IextM-IextM(end)), '-o')
xlabel('order')
title('change')

save_image("integr_extr")
%%
figure(1)
grid on

subplot(1,2,1)
hold on
plot(hN, I, '-o')
xlabel('h')
title('integral')

subplot(1,2,2)
loglog(hN, errN, '-o')
xlabel('h')
title('error')

%% Check order of convergence
figure(1)
hold on
loglog(hN, hN.^4)
legend('I_{err}', 'O(h^4)')
hold off

save_image("integr_conv")
close all