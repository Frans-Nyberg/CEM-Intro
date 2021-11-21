Nt = 8192; % Number of time steps
%Lx = .05; Ly = .04; Lz = .03; % Cavity dimensions in meters
%Nx = 25; Ny = 20; Nz = 15; % Number of cells along each axis
L = 0.05; % Cavity dimensions in meters
N = 37; % Number of cells along each axis

%% Run and sample
%[Et, Dt, c0] = Wave3D(Nt, Lx, Ly, Lz, Nx, Ny, Nz);
s1 = round(N*pi/10);
s2 = round(N*sqrt(3)/10);
s3 = round(N*7/11);
[Et, Dt, c0] = Wave3D(Nt, L, L, L, N, N, N, s1,s2,s3);

%%
t1 = Nt * Dt;
time = linspace(0,t1,Nt);
omega = 2*pi/t1*linspace(0, Nt-1, Nt);

%% Field
plot(time(end-50:end), Et(end-50:end,1))

%% columnwise fft
Nmax = 300;
spectr = fft(Et);
amplitude = abs(spectr);
plot(omega(2:Nmax), amplitude(2:Nmax,:))
xlabel("omega [rad/s]")
ylabel("amplitude")

%% Analytical reference
M = 3;  N = 4;  P = 4;
an_modes = 0.5*c0/L * sqrt(sum(make_trips(M,N,P).^2,2));
hold on
for k=1:length(an_modes)
    xline(2*pi*an_modes(k));
end

%% Save
savedir = 'results_2_2/';
save_image = @(filename) saveas(gca, strcat(savedir, filename), 'png');
save_image("cube_spectrum")