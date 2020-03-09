@everywhere using DelimitedFiles
@everywhere using SharedArrays
@everywhere include("constants.jl")

mutable struct Monitor
    totpop::SharedArray{Float64, 1}
    totnec::SharedArray{Float64, 1}
    vol::SharedArray{Float64, 1}
    Rvol::SharedArray{Float64, 1}
    totnew::SharedArray{Float64, 1}
    Rtotnew::SharedArray{Float64, 1}
    Shannon::SharedArray{Float64, 1}
    Simpson::SharedArray{Float64, 1}
    pops::SharedArray{Float64, 2}     # Total cell number per voxel (space)
    popt::SharedArray{Int64, 3}       # All populations cell number per time
    Vol2::SharedArray{Float64, 1}
    elapsed::Float64
    evalstep::Int64
    voxPop::Int64

    function Monitor(c::Constants)
        Neval = c.Neval
        N = c.N
        # Create monitor variables
        totpop = SharedArray{Float64, 1}((Neval,), init=S ->
                                                    setindex!(S, [c.P0], [1]))
        totnec = SharedArray{Float64, 1}((Neval,))
        vol = SharedArray{Float64, 1}((Neval,))
        Rvol = SharedArray{Float64, 1}((Neval,), init=S ->
                                                    setindex!(S, [1], [1]))
        totnew = SharedArray{Float64, 1}((Neval,))
        Rtotnew = SharedArray{Float64, 1}((Neval,))
        Shannon = SharedArray{Float64, 1}((Neval,))
        Simpson = SharedArray{Float64, 1}((Neval,), init=S ->
                                                    setindex!(S, [1], [1]))
        # Total cell number per voxel (space)
        function _init_pops!(S::SharedArray)
            S[1, 1] = c.P0
        end
        pops = SharedArray{Float64, 2}((2^c.alt, Neval), init=_init_pops!)
        # All populations cell number per time
        function _init_popt!(S::SharedArray)
            S[Int64(N / 2), Int64(N / 2), Int64(N / 2)] = c.P0
        end
        popt = SharedArray{Int64, 3}((N, N, N), init=_init_popt!)
        Vol2 = SharedArray{Float64, 1}((Neval,))
        # Parameters of evolving system
        elapsed = 0.0
        evalstep = 1
        voxPop = 0
        new(totpop, totnec, vol, Rvol, totnew, Rtotnew, Shannon, Simpson, pops,
            popt, Vol2, elapsed, evalstep, voxPop)
    end
end

function update_monitor_populations!(m::Monitor, c::Constants,
        G::SharedArray{Float64, 4}, Nec::SharedArray{Float64, 3},
        Act::SharedArray{Float64, 3}, i::Int64, j::Int64, k::Int64)
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
    m.Vol2[m.evalstep + 1] = size(findall(x -> x > c.threshold, m.popt), 1)
    for e = 1 : 2^c.alt
        if m.pops[e, m.evalstep + 1] > 0
            m.Shannon[m.evalstep + 1] = m.Shannon[m.evalstep + 1] -
            (m.pops[e, m.evalstep + 1] / m.totpop[m.evalstep + 1]) *
            log(m.pops[e, m.evalstep + 1] / m.totpop[m.evalstep + 1])
            m.Simpson[m.evalstep + 1] = m.Simpson[m.evalstep + 1] +
            (m.pops[e, m.evalstep + 1] / m.totpop[m.evalstep + 1])^2
        end
    end
    m.elapsed = time() - c.TimeStart
    m.evalstep = m.evalstep + 1
end

function monitor2files(m::Monitor, subfolder::String="files/")
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
