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
m = Monitor(c)

start = time()

open("files/Params.txt", "w") do file
    println(file, c.Grate, " ", c.Drate, " ", c.Mutrate, " ", c.Migrate)
end

# t = 0

for t in 1:c.Nstep
# @time while Vol2[evalstep] < 100000
    # t = t + 1;
    # Take care of local scope. Variables updated inside
    # for loop need to be assigned to global scope

    # global t;
    # global start;

    grid_time_step!(g, c, m, t)

    # Housekeeping
    if t % round(c.Nstep / c.Neval) == 0
        update_monitor_stats!(m, c)
        save_gen_space(g, t, c.N, "files/")
        print_curr_stats(m, t)

        global start = time()
    end
end

# Store tracking variables into files in `files` subfolder
monitor2files(m, "files/")
