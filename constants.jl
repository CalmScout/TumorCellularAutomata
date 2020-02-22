"""
    Contains constants / initial conditions of simulation.
"""
# value of random number generator seed
seedVal = 42;

# Set time of simulation
deltat = 4; # hours
# tspan = 6*10*365*deltat;
tspan = 1000;
Nstep = tspan/deltat;

# Set number of alterations
alt = 3;

# Create array of population cells, necrotics, and activity per voxel
N = 80;

G = zeros(N,N,N,2^alt);
Nec = zeros(N,N,N);
Act = zeros(N,N,N);
Rho = zeros(N,N,N);

# Set carrying capacity and initial cell number
P0 = 1e1;
K = 2e5;
# Assign initial cell number to population 1 at central voxel
G[Int64(N/2),Int64(N/2),Int64(N/2),1] = P0;

# Set number of time steps to be taken between system evaluations
Neval = Int64(ceil(Nstep / 20)) + 1;

# Create monitor variables
threshold = 0.2*K;
totpop = zeros(Neval);
totpop[1] = P0;
totnec = zeros(Neval);
totnec[1] = 0;
vol = zeros(Neval);
Rvol = zeros(Neval);
vol[1] = 0;
Rvol[1] = 1;
totnew = zeros(Neval);
Rtotnew = zeros(Neval);
totnew[1] = 0;
Rtotnew[1] = 0;
Shannon = zeros(Neval);
Simpson = zeros(Neval);
Shannon[1] = 0;
Simpson[1] = 1;
pops = zeros(2^alt,Neval); # Total cell number per voxel (space)
pops[1,1] = P0;
popt = Array{Int64}(undef,N,N,N); # All populations cell number per time
popt[Int64(N/2),Int64(N/2),Int64(N/2)] = P0;
Vol2 = zeros(Neval);
Vol2[1] = 0;
start = time();

# Create swapping matrix
Gnext = G;
G2 = G;
Necnext = Nec;
Actnext = Act;
Rhonext = Rho;

# Create swapping matrix
Gnext = G;
G2 = G;
Necnext = Nec;
Actnext = Act;
Rhonext = Rho;


MinGrate = 80;
MaxGrate = 250;    # Around 15 days
MinDrate = 80;
MaxDrate = 400;
MinMutrate = 80;
MaxMutrate = 240;
MinMigrate = 80;
MaxMigrate = 300;  # Around 25 days

Grate = 1;
Migrate = 10;
