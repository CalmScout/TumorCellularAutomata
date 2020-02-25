using DelimitedFiles

include("constants.jl")

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
    elapsed::Float64
    evalstep::Int64
    voxPop::Int64

    function Monitor(c::Constants)
        Neval = c.Neval
        N = c.N
        # Create monitor variables
        totpop = zeros(Neval)
        totpop[1] = c.P0
        totnec = zeros(Neval)
        vol = zeros(Neval)
        Rvol = zeros(Neval)
        Rvol[1] = 1
        totnew = zeros(Neval)
        Rtotnew = zeros(Neval)
        Shannon = zeros(Neval)
        Simpson = zeros(Neval)
        Simpson[1] = 1
        pops = zeros(2^c.alt, Neval) # Total cell number per voxel (space)
        pops[1, 1] = c.P0
        popt = Array{Int64}(undef, N, N, N) # All populations cell number per time
        popt[Int64(N / 2), Int64(N / 2), Int64(N / 2)] = c.P0
        Vol2 = zeros(Neval)
        # Parameters of evolving system
        elapsed = 0
        evalstep = 1
        voxPop = 0
        new(totpop, totnec, vol, Rvol, totnew, Rtotnew, Shannon, Simpson, pops,
            popt, Vol2, elapsed, evalstep, voxPop)
    end
end

function update_monitor_populations!(m::Monitor, c::Constants, G, Nec, Act,
    i::Int64, j::Int64, k::Int64)
    m.totpop[m.evalstep + 1] = m.totpop[m.evalstep + 1] + sum(G[i, j, k, :])
    m.totnec[m.evalstep + 1] = m.totnec[m.evalstep + 1] + sum(Nec[i, j, k])
    m.Rtotnew[m.evalstep + 1] = m.Rtotnew[m.evalstep + 1] + sum(Act[i, j, k])
    m.Rvol[m.evalstep + 1] = m.Rvol[m.evalstep + 1] + 1

    for e = 1 : 2^c.alt
        m.pops[e, m.evalstep + 1] = m.pops[e, m.evalstep + 1] + G[i, j, k, e]
    end

    if sum(G[i, j, k, :]) > c.threshold
        m.totnew[m.evalstep + 1] = m.totnew[m.evalstep + 1] + sum(Act[i, j, k])
        m.vol[m.evalstep + 1] = m.vol[m.evalstep + 1] + 1
    end
end

function update_monitor_stats!(m::Monitor, c::Constants)
    m.Shannon[m.evalstep + 1] = 0
    m.Simpson[m.evalstep + 1] = 0
    ROcc =  findall(x -> x > c.threshold, m.popt)
    m.Vol2[m.evalstep + 1] = size(ROcc, 1)
    ROcc = []
    for e = 1 : 2^c.alt
        if m.pops[e, m.evalstep + 1] > 0
            m.Shannon[m.evalstep + 1] = m.Shannon[m.evalstep + 1] -
            (m.pops[e, m.evalstep + 1] / m.totpop[m.evalstep + 1]) *
            log(m.pops[e, m.evalstep + 1] / m.totpop[m.evalstep + 1])
            m.Simpson[m.evalstep + 1] = m.Simpson[m.evalstep + 1] +
            (m.pops[e, m.evalstep + 1] / m.totpop[m.evalstep + 1])^2
        end
    end
    m.elapsed = m.elapsed + time() - start
    m.evalstep = m.evalstep + 1
end

function monitor2files(m::Monitor, subfolder="files/")
    dir_to_save = joinpath(@__DIR__, subfolder)
    writedlm(joinpath(dir_to_save, "Totpop.txt"), m.totpop)
    writedlm(joinpath(dir_to_save, "Totnec.txt"), m.totnec)
    writedlm(joinpath(dir_to_save, "VolPET.txt"), m.vol)
    writedlm(joinpath(dir_to_save, "Vol_real.txt"), m.Rvol)
    writedlm(joinpath(dir_to_save, "ActPET.txt"), m.totnew)
    writedlm(joinpath(dir_to_save, "Act_real.txt"), m.Rtotnew)
    writedlm(joinpath(dir_to_save, "Shannon.txt"), m.Shannon)
    writedlm(joinpath(dir_to_save, "Simpson.txt"), m.Simpson)
    writedlm(joinpath(dir_to_save, "Genspop.txt"), m.pops)
end
