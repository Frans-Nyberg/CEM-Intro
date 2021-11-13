%% following Ch 3 in T.Rylander 2013
savedir = 'results_2_2/';
save_image = @(filename) saveas(gca, strcat(savedir, filename), 'png');
a = pi;
nx = 30;
dx         = a/nx;           % The cell size 
d2tmax     = 1.9*dx;         % The time step must satisfy
                             % 2*dt < 2*dx for stability

%% Eigenfrequency study
time = 30*pi;
[omega, s1, s2] = Wave1D(a, time, nx, d2tmax);
plots(omega, s1,s2, "t=30\pi")
save_image("spect-t30")

time = 31*pi;
[omega, s1, s2] = Wave1D(a, time, nx, d2tmax);
plots(omega, s1,s2, "t=31\pi")
save_image("spect-t31")

close all

%% Stability study

time = 30*pi;
d2tmax = 2.1*dx;
[omega, s1, s2] = Wave1D(a, time, nx, d2tmax);
plots(omega, s1,s2, "t=30\pi")
save_image("stab-t30")

%%
function plots(omega,s1,s2,titl)
figure()
subplot(1,2,1)
plot(omega, abs(s1), '-x')
xlabel("\omega [rad/s]")
ylabel("|s_{a/5}|")
title(titl)
subplot(1,2,2)
plot(omega, abs(s2), '-x')
xlabel("\omega [rad/s]")
ylabel("|s_{a/2}|")
title(titl)
end