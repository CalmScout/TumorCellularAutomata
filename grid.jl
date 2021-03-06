include("constants.jl")

mutable struct Grid
    G::Array{Float64, 4}
    Nec::Array{Float64, 3}
    Act::Array{Float64, 3}
    Rho::Array{Float64, 3}
    Gnext::Array{Float64, 4}
    G2::Array{Float64, 4}
    Necnext::Array{Float64, 3}
    Actnext::Array{Float64, 3}
    Rhonext::Array{Float64, 3}
    Occ::Array{CartesianIndex{3}, 1}
    ROcc::Array{CartesianIndex{3}, 1}

    function Grid(c::Constants)
        # Create array of population cells, necrotics, and activity per voxel
        G = zeros(c.N, c.N, c.N, 2^c.alt)
        Nec = zeros(c.N, c.N, c.N)
        Act = zeros(c.N, c.N, c.N)
        Rho = zeros(c.N, c.N, c.N)
        # Assign initial cell number to population 1 at central voxel
        G[Int64(c.N / 2), Int64(c.N / 2), Int64(c.N / 2), 1] = c.P0
        # Create swapping matrices
        Gnext = copy(G)
        G2 = copy(G)
        Necnext = copy(Nec)
        Actnext = copy(Act)
        Rhonext = copy(Rho)
        Occ = [CartesianIndex(Int(c.N / 2), Int(c.N / 2), Int(c.N / 2))]
        ROcc = [CartesianIndex(Int(c.N / 2), Int(c.N / 2), Int(c.N / 2))]
        new(G, Nec, Act, Rho, Gnext, G2, Necnext, Actnext, Rhonext, Occ, ROcc)
    end
end

function grid_time_step!(g::Grid, c::Constants, m::Monitor, t::Int64)
    for l in 1:length(g.Occ)
        i = Int(g.Occ[l][1])
        j = Int(g.Occ[l][2])
        k = Int(g.Occ[l][3])

        # Reinitialize activity at each time step
        g.Act[i, j, k] = 0
        # Only evaluate voxel if there is at least 1 cell
        if sum(g.G[i, j, k, :]) > 0
            for e in 1:2^c.alt
                # Only evaluate population if there is at least 1 cell
                if g.G[i, j, k, e] > 0
                    # receive binary representation
                    binGb = decimal2binstr(e, c.alt)

                    # Retrieve voxel info
                    Popgen = g.G[i, j, k, e]
                    Popvox = sum(g.G[i, j, k, :])
                    Necvox = g.Nec[i, j, k]

                    # Reproduction event
                    born =
                    reproduction_event!(g, c, Popgen, Popvox, Necvox, i, j, k,
                    e, binGb)

                    # Death event
                    dead = death_event!(g, c, binGb, Popvox, Necvox, Popgen, i,
                    j, k, e)

                    # Migration event
                    migration_event!(g, c, binGb, Popvox, Popgen, Necvox, i, j,
                    k, e)

                    # Mutation event
                    mutation_event!(g, c, binGb, Popgen, i, j, k, e)

                end
            end

            if t % c.NstepNevalRatio == 0
                update_monitor_populations!(m, c, g.G, g.Nec, g.Act, i,
                    j, k)
            end
        end
    end

    #update swapping matrices
    g.G = copy(g.Gnext)
    g.Nec = copy(g.Necnext)
    g.Act = copy(g.Actnext)
    g.Rho = copy(g.Rhonext)

    g.G2 = sum(g.G, dims = 4)
    m.popt = g.G2[:, :, :, 1] + g.Nec
    g.Occ = findall(x -> x > 0, m.popt)
end

function normalize_prob(Prep::Float64)
    Prep = max(min(Prep, 1), 0)
end

function reproduction_event!(g::Grid, c::Constants, Popgen::Float64,
    Popvox::Float64, Necvox::Float64, i::Int64, j::Int64, k::Int64, e::Int64,
    binGb::Array{Float64, 1})
    """
        Reproduction event
    """
    grate = c.Grate * (1 - binGb' * c.Gweight)
    Prep = c.deltat / grate *(1-(Popvox + Necvox) / c.K)
    Prep = normalize_prob(Prep)
    born = rand(Binomial(Int64(Popgen), Prep))
    g.Gnext[i, j, k, e] = g.Gnext[i, j, k, e] + born
    g.Actnext[i, j, k] = g.Actnext[i, j, k] + born
    return born
end

function death_event!(g::Grid, c::Constants, binGb::Array{Float64, 1},
    Popvox::Float64, Necvox::Float64, Popgen::Float64, i::Int64, j::Int64,
    k::Int64, e::Int64)
    """
        Death event
    """
    drate = c.Drate * (1 - binGb' * c.Dweight)
    Pkill = c.deltat / drate * (Popvox + Necvox) / c.K
    Pkill = normalize_prob(Pkill)
    dead = rand(Binomial(Int64(Popgen), Pkill))
    g.Gnext[i, j, k, e] = g.Gnext[i, j, k, e] - dead
    g.Necnext[i, j, k] = g.Necnext[i, j, k] + dead
    return dead
end

function migration_event!(g::Grid, c::Constants, binGb::Array{Float64, 1},
    Popvox::Float64, Popgen::Float64, Necvox::Float64, i::Int64, j::Int64,
    k::Int64, e::Int64)
    """
        Migration event
    """
    migrate = c.Migrate * (1 - binGb' * c.Migweight)
    Pmig = c.deltat / migrate * (Popvox + Necvox) / c.K
    Pmig = normalize_prob(Pmig)
    migrants = rand(Binomial(Int64(Popgen), Pmig))
    neigh = 0
    moore = 26
    vonN = 6
    multinom = Multinomial(migrants, c.wcube)
    gone = rand(multinom)
    N = c.N
    for movi in [-1, 0, 1]
        for movj in [-1, 0, 1]
            for movk = [-1, 0, 1]
                xmov = i + movi
                ymov = j + movj
                zmov = k + movk
                if xmov < N+1 && ymov < N+1 && zmov < N+1 && xmov > 0 &&
                    ymov > 0 && zmov > 0 && abs(movi)+abs(movj)+abs(movk)!=0
                    neigh = neigh + 1
                    g.Gnext[xmov, ymov, zmov, e] =
                                g.Gnext[xmov, ymov, zmov, e] + gone[neigh]
                    g.Gnext[i, j, k, e] = g.Gnext[i, j, k, e] - gone[neigh]
                end
            end
        end
    end
end

function mutation_event!(g::Grid, c::Constants, binGb::Array{Float64, 1},
    Popgen::Float64, i::Int64, j::Int64, k::Int64, e::Int64)
    """
        Mutation event
    """
    mutrate = c.Mutrate * (1 - binGb' * c.Mutweight)
    Pmut = c.deltat / mutrate * (Popgen / c.K)
    Pmut = normalize_prob(Pmut)
    r = rand(1)
    r = r[1]
    if r < Pmut && e != 2^c.alt
        # Pick a random empty slot and turn it to mutated
        nonalter = findall(x -> x < 1, binGb)
        r2 = rand(1:length(nonalter))
        mutating = nonalter[r2]
        binGb[mutating] = 1

        # Switch binary array back to binary string
        binGc = string(Int(binGb[1]), Int(binGb[2]), Int(binGb[3]))

        # Code below retrieves back decimal number from binary string
        decG = parse(Int, binGc, base=2) + 1
        g.Gnext[i, j, k, e] = g.Gnext[i, j, k, e] - 1
        g.Gnext[i, j, k, decG] = g.Gnext[i, j, k, decG] + 1
    end
end
