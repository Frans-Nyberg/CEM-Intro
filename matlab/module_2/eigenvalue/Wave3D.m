% From book of Anders Bondeson, Thomas Rylander, Pär Ingelström
% in the 2005 version
function [Et, Dt, c0] = Wave3D(Nt, Lx, Ly, Lz, Nx, Ny, Nz, s1, s2, s3)
% Physical constants
eps0 = 8.8541878e-12; % Permittivity of vacuum
mu0 = 4e-7 * pi; % Permeability of vacuum
c0 = 299792458; % Speed of light in vacuum
% Parameter initiation
Cx = Nx / Lx; % Inverse cell dimensions
Cy = Ny / Ly; 
Cz = Nz / Lz; 
Dt = 1/(c0*norm([Cx Cy Cz])); % Time step
% Allocate field matrices
Ex = zeros(Nx , Ny+1, Nz+1);
Ey = zeros(Nx+1, Ny , Nz+1);
Ez = zeros(Nx+1, Ny+1, Nz );
Hx = zeros(Nx+1, Ny , Nz );
Hy = zeros(Nx , Ny+1, Nz );
Hz = zeros(Nx , Ny , Nz+1);
% Allocate time signals
Et = zeros(Nt,3);
% Initiate fields with noise (except on the boundary)
Ex( : , 2:Ny, 2:Nz) = rand(Nx , Ny-1, Nz-1) - 0.5;
Ey(2:Nx, : , 2:Nz) = rand(Nx-1, Ny , Nz-1) - 0.5;
Ez(2:Nx, 2:Ny, : ) = rand(Nx-1, Ny-1, Nz ) - 0.5;
% Time stepping
for n = 1:Nt
% Update H everywhere
Hx = Hx + (Dt/mu0)*(diff(Ey,1,3)*Cz - diff(Ez,1,2)*Cy);
Hy = Hy + (Dt/mu0)*(diff(Ez,1,1)*Cx - diff(Ex,1,3)*Cz);
Hz = Hz + (Dt/mu0)*(diff(Ex,1,2)*Cy - diff(Ey,1,1)*Cx);
% Update E everywhere except on boundary
Ex(:,2:Ny,2:Nz) = Ex(:,2:Ny,2:Nz) + (Dt /eps0) * ...
(diff(Hz(:,:,2:Nz),1,2)*Cy - diff(Hy(:,2:Ny,:),1,3)*Cz);
Ey(2:Nx,:,2:Nz) = Ey(2:Nx,:,2:Nz) + (Dt /eps0) * ...
(diff(Hx(2:Nx,:,:),1,3)*Cz - diff(Hz(:,:,2:Nz),1,1)*Cx);
Ez(2:Nx,2:Ny,:) = Ez(2:Nx,2:Ny,:) + (Dt /eps0) * ...
(diff(Hy(:,2:Ny,:),1,1)*Cx - diff(Hx(2:Nx,:,:),1,2)*Cy);
% Sample the electric field at chosen points
Et(n,:) = [Ex(s1,s2,s3) Ey(s1,s2,s3) Ez(s1,s2,s3)];
end