# Will come in handy to random sample from binomial distributions
using Distributions
using Random
using DelimitedFiles

# import predefined constants / parameters
include("constants.jl")
# import helper functions
include("tools.jl")

Random.seed!(seedVal)

# Create array of population cells, necrotics, and activity per voxel
G = zeros(N, N, N, 2^alt)
Nec = zeros(N, N, N)
Act = zeros(N, N, N)
Rho = zeros(N, N, N)

# Assign initial cell number to population 1 at central voxel
G[Int64(N/2), Int64(N/2), Int64(N/2), 1] = P0

# Create swapping matrix
Gnext = G
G2 = G
Necnext = Nec
Actnext = Act
Rhonext = Rho

# Create monitor variables
totpop = zeros(Neval)
totpop[1] = P0
totnec = zeros(Neval)
vol = zeros(Neval)
Rvol = zeros(Neval)
Rvol[1] = 1
totnew = zeros(Neval)
Rtotnew = zeros(Neval)
Shannon = zeros(Neval)
Simpson = zeros(Neval)
Simpson[1] = 1
pops = zeros(2^alt, Neval) # Total cell number per voxel (space)
pops[1, 1] = P0
popt = Array{Int64}(undef, N, N, N) # All populations cell number per time
popt[Int64(N/2), Int64(N/2), Int64(N/2)] = P0
Vol2 = zeros(Neval)
start = time()

GrateInit = 1
MigrateInit = 10
GrateInit, MigrateInit = adjust_grate_migrate(GrateInit, MigrateInit,
                                MinGrate, MaxGrate, MinMigrate, MaxMigrate)
const Grate = GrateInit
const Migrate = MigrateInit
const Drate = rand(Uniform(MinDrate, MaxDrate))
const Mutrate = rand(Uniform(MinMutrate, MaxMutrate))

open("files/Params.txt", "w") do file
    println(file, Grate, " ", Drate, " ", Mutrate, " ", Migrate)
end

# Create weights for surrounding voxels (Moore neighbourhood)
c = 0
wcube = zeros(26)
sumcube = 0

for i in [-1, 0, 1]
    for j in [-1, 0, 1]
        for k in [-1, 0, 1]
            global c
            global wcube
            global sumcube
            if abs(i) + abs(j) + abs(k) != 0
                c = c + 1
                wcube[c] = 1 / sqrt(abs(i) + abs(j) + abs(k))
                sumcube = sumcube + wcube[c]
            end
        end
    end
end

c = 0
for i in [-1, 0, 1]
    for j in [-1, 0, 1]
        for k in [-1, 0, 1]
            global c
            global wcube
            global sumcube
            if abs(i) + abs(j) + abs(k) != 0
                c = c + 1
                wcube[c] = wcube[c] / sumcube
            end
        end
    end
end

# Let the system evolve
elapsed = 0
evalstep = 1
voxPop = 0
Occ = [CartesianIndex(Int(N/2), Int(N/2), Int(N/2))]
ROcc = [CartesianIndex(Int(N/2), Int(N/2), Int(N/2))]
# t = 0

for t in 1:Nstep
# @time while Vol2[evalstep] < 100000
    # t = t + 1;
    # Take care of local scope. Variables updated inside
    # for loop need to be assigned to global scope
    global G
    # global t;
    global Nec
    global Act
    global Rho
    global Gnext
    global Necnext
    global Actnext
    global Rhonext
    global evalstep
    global elapsed
    global totpop
    global totnec
    global vol
    global Rvol
    global totnew
    global Rtotnew
    global Shannon
    global Simpson
    global pops
    global popt
    global Occ
    global G2
    global Gweight
    global Dweight
    global Migweight
    global Mutweight
    global wcube
    global Vol2
    global ROcc
    # global start;

    for l in 1:length(Occ)
        i = Int(Occ[l][1])
        j = Int(Occ[l][2])
        k = Int(Occ[l][3])

        # Reinitialize activity at each time step
        Act[i,j,k] = 0
        # Only evaluate voxel if there is at least 1 cell
        if sum(G[i, j, k, :]) > 0
            for e in 1:2^alt
                # Only evaluate population if there is at least 1 cell
                if G[i, j, k, e] > 0
                    # receive binary representation
                    binGb = decimal2binstr(e)

                    # Retrieve voxel info
                    Popgen = G[i, j, k, e]
                    Popvox = sum(G[i, j, k, :])
                    Necvox = Nec[i, j, k]

                    # Reproduction event
                    born =
                    reproduction_event(Popgen, Popvox, K, Necvox, Grate,
                                                deltat, i, j, k, e, alt, binGb)

                    # Death event
                    dead = death_event(Drate, binGb, Dweight, deltat,
                                Popvox, Necvox, K, Popgen, i, j, k, e)

                    # Migration event
                    migration_event(Migrate, binGb, Migweight, deltat,
                                    Popvox, Popgen, Necvox, K, i, j, k, e)

                    # Mutation event
                    mutation_event(Mutrate, binGb, Mutweight, deltat, Popgen, K)
                end
            end

            # Housekeeping
            if t % round(Nstep / Neval) == 0
                totpop[evalstep + 1] = totpop[evalstep + 1] + sum(G[i, j, k, :])
                totnec[evalstep + 1] = totnec[evalstep + 1] + sum(Nec[i, j, k])
                Rtotnew[evalstep + 1] = Rtotnew[evalstep + 1] + sum(Act[i, j, k])
                Rvol[evalstep + 1] = Rvol[evalstep + 1] + 1

                for e = 1 : 2^alt
                    pops[e, evalstep + 1] = pops[e, evalstep + 1] + G[i, j, k, e]
                end

                if sum(G[i, j, k, :]) > threshold
                    totnew[evalstep + 1] = totnew[evalstep + 1] + sum(Act[i, j, k])
                    vol[evalstep + 1] = vol[evalstep + 1] + 1
                end
            end
        end
    end

    G = Gnext
    Nec = Necnext
    Act = Actnext
    Rho = Rhonext
    G2 = sum(G, dims = 4)
    popt = G2[:, :, :, 1] + Nec

    Occ = findall(x -> x > 0, popt)
    # Housekeeping
    if t % round(Nstep / Neval) == 0
        Shannon[evalstep + 1] = 0
        Simpson[evalstep + 1] = 0
        ROcc =  findall(x -> x > threshold, popt)
        Vol2[evalstep + 1] = size(ROcc, 1)
        ROcc = []
        for e = 1 : 2^alt
            if pops[e, evalstep + 1] > 0
                Shannon[evalstep + 1] = Shannon[evalstep + 1] -
                (pops[e, evalstep + 1] / totpop[evalstep + 1]) *
                log(pops[e, evalstep + 1] / totpop[evalstep + 1])
                Simpson[evalstep + 1] = Simpson[evalstep + 1] +
                (pops[e, evalstep + 1] / totpop[evalstep + 1])^2
            end
        end

        save_gen_space(G, Act, Nec, t, N, "files/")

        elapsed = elapsed + time() - start
        println("Cell no: ", totpop[evalstep + 1], "; Volume: ", Rvol[evalstep +
        1], "; Activity: ", Rtotnew[evalstep + 1], "; Necrotics: ",
        totnec[evalstep + 1], "; Het: ", Shannon[evalstep + 1])
        println("     Volume (alt): ", Vol2[evalstep + 1])
        println("Time step: ", t, "; Time elapsed: ", elapsed)
        evalstep = evalstep + 1
        global start = time()
    end
end

# Store tracking variables into files in `files` subfolder
dir_to_save = joinpath(@__DIR__, "files/")
writedlm(joinpath(dir_to_save, "Totpop.txt"), totpop)
writedlm(joinpath(dir_to_save, "Totnec.txt"), totnec)
writedlm(joinpath(dir_to_save, "VolPET.txt"), vol)
writedlm(joinpath(dir_to_save, "Vol_real.txt"), Rvol)
writedlm(joinpath(dir_to_save, "ActPET.txt"), totnew)
writedlm(joinpath(dir_to_save, "Act_real.txt"), Rtotnew)
writedlm(joinpath(dir_to_save, "Shannon.txt"), Shannon)
writedlm(joinpath(dir_to_save, "Simpson.txt"), Simpson)
writedlm(joinpath(dir_to_save, "Genspop.txt"), pops)
