mutable struct Monitor
    totpop::Array{Float64, 1}
    totnec::Array{Float64, 1}
    vol::Array{Float64, 1}
    Rvol::Array{Float64, 1}
    totnew::Array{Float64, 1}
    Rtotnew::Array{Float64, 1}
    Shannon::Array{Float64, 1}
    Simpson::Array{Float64, 1}
    pops::Array{Float64, 2}     # Total cell number per voxel (space)
    popt::Array{Int64, 3}       # All populations cell number per time
    Vol2::Array{Float64, 1}
    function Monitor(Neval::Int64)
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
        new(totpop, totnec, vol, Rvol, totnew, Rtotnew, Shannon, Simpson, pops,
            popt, Vol2)
    end
end

function update_monitor_populations!(m::Monitor, evalstep::Int64,
    G, Nec, Act, i::Int64, j::Int64, k::Int64, threshold::Float64)
    m.totpop[evalstep + 1] = m.totpop[evalstep + 1] + sum(G[i, j, k, :])
    m.totnec[evalstep + 1] = m.totnec[evalstep + 1] + sum(Nec[i, j, k])
    m.Rtotnew[evalstep + 1] = m.Rtotnew[evalstep + 1] + sum(Act[i, j, k])
    m.Rvol[evalstep + 1] = m.Rvol[evalstep + 1] + 1

    for e = 1 : 2^alt
        m.pops[e, evalstep + 1] = m.pops[e, evalstep + 1] + G[i, j, k, e]
    end

    if sum(G[i, j, k, :]) > threshold
        m.totnew[evalstep + 1] = m.totnew[evalstep + 1] + sum(Act[i, j, k])
        m.vol[evalstep + 1] = m.vol[evalstep + 1] + 1
    end
end

function update_monitor_stats!(m::Monitor, evalstep::Int64,
    threshold::Float64)
    m.Shannon[evalstep + 1] = 0
    m.Simpson[evalstep + 1] = 0
    ROcc =  findall(x -> x > threshold, m.popt)
    m.Vol2[evalstep + 1] = size(ROcc, 1)
    ROcc = []
    for e = 1 : 2^alt
        if m.pops[e, evalstep + 1] > 0
            m.Shannon[evalstep + 1] = m.Shannon[evalstep + 1] -
            (m.pops[e, evalstep + 1] / m.totpop[evalstep + 1]) *
            log(m.pops[e, evalstep + 1] / m.totpop[evalstep + 1])
            m.Simpson[evalstep + 1] = m.Simpson[evalstep + 1] +
            (m.pops[e, evalstep + 1] / m.totpop[evalstep + 1])^2
        end
    end
end
