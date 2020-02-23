"""
    Contains constants / initial conditions of simulation.
"""
# value of random number generator seed
const seedVal = 42;

# Set time of simulation
const deltat = 4; # hours
# tspan = 6*10*365*deltat;
const tspan = 1000;
const Nstep = tspan/deltat;

# Set number of alterations
const alt = 3;

# Grid dimensions N x N x N
const N = 80;

# Set carrying capacity and initial cell number
const P0 = 1e1;
const K = 2e5;

# Set number of time steps to be taken between system evaluations
const Neval = Int64(ceil(Nstep / 20)) + 1;

const threshold = 0.2 * K;

const MinGrate = 80;
const MaxGrate = 250;    # Around 15 days
const MinDrate = 80;
const MaxDrate = 400;
const MinMutrate = 80;
const MaxMutrate = 240;
const MinMigrate = 80;
const MaxMigrate = 300;  # Around 25 days

# Set all weights (slightly tuned)
const Gweight = [0.32, 0.28, 0.25];
const Dweight = [-0.15, -0.05, -0.45];
const Mutweight = [0.18, 0.18, 0.32];
const Migweight = [0.65, 0.05, 0.05];
