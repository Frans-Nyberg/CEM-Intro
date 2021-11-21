Nt = 8192; % Number of time steps
%Lx = .05; Ly = .04; Lz = .03; % Cavity dimensions in meters
%Nx = 25; Ny = 20; Nz = 15; % Number of cells along each axis
L = 0.05; % Cavity dimensions in meters
Nh = 37; % Number of cells along each axis
c0 = 299792458; % Speed of light in vacuum

%% Run and sample
%[Et, Dt, c0] = Wave3D(Nt, Lx, Ly, Lz, Nx, Ny, Nz);
r1 = 1/7; r2 = sqrt(5)/5; r3 = 8/11;
s1 = round(Nh*r1);
s2 = round(Nh*r2);
s3 = round(Nh*r3);
dt = L/Nh/(c0*sqrt(3));
%%
[Et, Dt] = Wave3D(Nt, L, L, L, Nh, Nh, Nh, s1,s2,s3,dt);

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
legend("Ex","Ey","Ez", 'Autoupdate', 'off')
side = linspace(-1/2, 1/2, Nh);
title(sprintf("at (%.4f, %.4f, %.4f)L", side(s1), side(s2), side(s3)))

%% Analytical reference
M = 4;  N = 4;  P = 5;
an_modes = 0.5*c0/L * sqrt(sum(make_trips(M,N,P).^2,2));
hold on
for k=1:length(an_modes)
    if 2*pi*an_modes(k) <= omega(Nmax)
        xline(2*pi*an_modes(k));
    end
end

%% Save
savedir = 'results_2_2/';
save_image = @(filename) saveas(gca, strcat(savedir, filename), 'png');
save_image("cube_spectrum")