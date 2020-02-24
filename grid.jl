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

    function Grid(N::Int64 , alt::Int64, P0::Float64)
        # Create array of population cells, necrotics, and activity per voxel
        G = zeros(N, N, N, 2^alt)
        Nec = zeros(N, N, N)
        Act = zeros(N, N, N)
        Rho = zeros(N, N, N)
        # Assign initial cell number to population 1 at central voxel
        G[Int64(N/2), Int64(N/2), Int64(N/2), 1] = P0
        # Create swapping matrices
        Gnext = G
        G2 = G
        Necnext = Nec
        Actnext = Act
        Rhonext = Rho
        new(G, Nec, Act, Rho, Gnext, G2, Necnext, Actnext, Rhonext)
    end
end

function grid_time_step!(g::Grid, m::Monitor, t, Occ, alt, K, Grate, Drate, Dweight,
    Migweight, Mutrate, Mutweight, deltat)
    for l in 1:length(Occ)
        i = Int(Occ[l][1])
        j = Int(Occ[l][2])
        k = Int(Occ[l][3])

        # Reinitialize activity at each time step
        g.Act[i, j, k] = 0
        # Only evaluate voxel if there is at least 1 cell
        if sum(g.G[i, j, k, :]) > 0
            for e in 1:2^alt
                # Only evaluate population if there is at least 1 cell
                if g.G[i, j, k, e] > 0
                    # receive binary representation
                    binGb = decimal2binstr(e)

                    # Retrieve voxel info
                    Popgen = g.G[i, j, k, e]
                    Popvox = sum(g.G[i, j, k, :])
                    Necvox = g.Nec[i, j, k]

                    # Reproduction event
                    born =
                    reproduction_event!(g, Popgen, Popvox, K, Necvox, Grate,
                                                deltat, i, j, k, e, alt, binGb)

                    # Death event
                    dead = death_event!(g, Drate, binGb, Dweight, deltat,
                                Popvox, Necvox, K, Popgen, i, j, k, e)

                    # Migration event
                    migration_event!(g, Migrate, binGb, Migweight, deltat,
                                    Popvox, Popgen, Necvox, K, i, j, k, e)

                    # Mutation event
                    mutation_event!(g, Mutrate, binGb, Mutweight, deltat, Popgen, K)
                end
            end

            # Housekeeping
            if t % round(Nstep / Neval) == 0
                update_monitor_populations!(m, evalstep, g.G, g.Nec, g.Act, i,
                    j, k, threshold)
            end
        end
    end
end
