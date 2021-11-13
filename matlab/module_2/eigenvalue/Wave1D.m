% --------------------------------------------------------------
% Time step 1D wave equation using two time-levels f0 & f1
% --------------------------------------------------------------
function [omega, s1, s2] = Wave1D(a, time, nx, d2tmax)

% Arguments:
%    a     = the length of the interval
%    time  = the total time interval for the simulation
%    nx    = the number of subintervals in the domain (0,a)
% Returns:
%    omega = the angular frequencies
%    s1    = the complex Fourier transform of data at x = a/5
%    s2    = the complex Fourier transform of data at x = a/2

f0         = randn(nx+1, 1); % Initialize with random numbers
f0(1,1)    = 0;              % Boundary condition at x = 0
f0(nx+1,1) = 0;              % Boundary condition at x = a

f1         = randn(nx+1, 1); % Initialize with random numbers
f1(1,1)    = 0;              % Boundary condition at x = 0
f1(nx+1,1) = 0;              % Boundary condition at x = a

dx         = a/nx;           % The cell size 

ntime = round(time/d2tmax + 1);  % The number of time steps
dt = time/(2*ntime);             % The time step

% Initialize the coefficient matrix for updating the solution f
A = spalloc(nx+1,nx+1,3*(nx+1)); % Sparse empty matrix with 
                                 % 3*(nx+1) nonzero entries
for i = 2:nx
  A(i,i)   = 2*(1-(dt/dx)^2);    % Diagonal entries
  A(i,i+1) = (dt/dx)^2;          % Upper diagonal entries
  A(i,i-1) = (dt/dx)^2;          % Lower diagonal entries
end

% Time step and sample the solution
% Sample location #1 is close to the left boundary
% Sample location #2 is at the midpoint of the domain
for itime = 1:ntime % Every 'itime' means two time steps 'dt'

  f0               = A*f1 - f0;          % Update
  sign1(2*itime-1) = f0(round(1+nx/5));  % Sample at location #1
  sign2(2*itime-1) = f0(round(1+nx/2));  % Sample at location #2

  f1               = A*f0 - f1;          % Update               
  sign1(2*itime)   = f1(round(1+nx/5));  % Sample at location #1
  sign2(2*itime)   = f1(round(1+nx/2));  % Sample at location #2

end

% Compute the discrete Fourier transform of 
% the time-domain signals
spectr1     = fft(sign1); 
spectr2     = fft(sign2);

% In the MATLAB implementation of the function fft(), 
% the first half of the output corresponds to positive frequency
s1(1:ntime) = spectr1(1:ntime); 
s2(1:ntime) = spectr2(1:ntime);

% Frequency vector for use with 's1' and 's2'
omega       = (2*pi/time)*linspace(0, ntime-1, ntime);

