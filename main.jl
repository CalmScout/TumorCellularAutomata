# Will come in handy to random sample from binomial distributions
using Distributions
using Random

# import predefined constants / parameters
include("constants.jl")
# import helper functions
include("tools.jl")
# main data structure
include("grid.jl")

const seedVal = 42
Random.seed!(seedVal)

c = Constants()
g = Grid(c)

# Create monitor variable
m = Monitor(c)

start = time()

open("files/Params.txt", "w") do file
    println(file, c.Grate, " ", c.Drate, " ", c.Mutrate, " ", c.Migrate)
end

# Let the system evolve
elapsed = 0
evalstep = 1
voxPop = 0
Occ = [CartesianIndex(Int(c.N / 2), Int(c.N / 2), Int(c.N / 2))]
ROcc = [CartesianIndex(Int(c.N / 2), Int(c.N / 2), Int(c.N / 2))]
# t = 0

for t in 1:c.Nstep
# @time while Vol2[evalstep] < 100000
    # t = t + 1;
    # Take care of local scope. Variables updated inside
    # for loop need to be assigned to global scope

    # global t;

    global evalstep
    global elapsed

    global Occ
    global Gweight
    global Dweight
    global Migweight
    global Mutweight
    global wcube
    global ROcc
    # global start;

    grid_time_step!(g, c, m, t, Occ)

    m.popt = g.G2[:, :, :, 1] + g.Nec

    Occ = findall(x -> x > 0, m.popt)
    # Housekeeping
    if t % round(c.Nstep / c.Neval) == 0
        update_monitor_stats!(m, c, evalstep)
        save_gen_space(g, t, c.N, "files/")

        elapsed = elapsed + time() - start
        print_curr_stats(m, t, elapsed, evalstep)
        evalstep = evalstep + 1

        global start = time()
    end
end

# Store tracking variables into files in `files` subfolder
monitor2files(m, "files/")
