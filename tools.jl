"""
    Module contains tools for testing, saving results to files, visualizastion
    etc.
"""

include("monitor.jl")
include("grid.jl")

function compare_files(path1, path2)
    """
    Reads two text files from provided pathes and compares them.
    """
    return read(path1, String) == read(path2, String)
end

function compare_folders(path_folder_1, path_folder_2)
    """
    Compare all files in folders `path_folder_1` and `path_folder_2`, returns
    true if number and content of files are equal.
    """
    files_path_1 = readdir(path_folder_1)
    files_path_2 = readdir(path_folder_2)
    # check id we have the same number of files in both of the folders
    if length(files_path_1) != length(files_path_2)
        return false
    end
    # try catches error if we do not have file in the second folder but have
    # in the first one
    try
        for filename in files_path_1
            file_1_path = joinpath(path_folder_1, filename)
            file_2_path = joinpath(path_folder_2, filename)
            # check if files with the same names have the same content
            curr_file_comparison = compare_files(file_1_path, file_2_path)
            if !curr_file_comparison
                return false
            end
        end
    catch
        e
        if isa(e, LoadError)
            return false
        end
    end
    return true
end

function decimal2binstr(e)
    """
        Converts decimal number to binary string.
    """

    binG = string(e-1, base=2);     # First of all, we retrieve binary string
    while length.(binG) < 3         # As the string may not be of length 0, we solve this
        binG = string("0",binG);
    end
    binGa = split(binG,"");         # Now we need an array instead of a string
    binGb = zeros(length(binGa));   # We create a new variable to store int array
    for i in 1:length(binGa)        # We convert array elements from char to int
        binGb[i] = parse(Int,binGa[i]);
    end
    return binGb
end

function reproduction_event!(g, Popgen, Popvox, K, Necvox, Grate, deltat, i, j, k,
    e, alt, binGb)
    """
        Reproduction event
    """
    grate = Grate*(1-binGb[1]*Gweight[1]-binGb[2]*Gweight[2]-binGb[3]*Gweight[3]);
    Prep = (deltat/grate)*(1-(Popvox+Necvox)/K);
    if Prep > 1
        Prep = 1;
    end
    if Prep < 0
        Prep = 0;
    end
    born = rand(Binomial(Int64(Popgen), Prep));
    g.Gnext[i,j,k,e] = g.Gnext[i,j,k,e] + born;
    g.Actnext[i,j,k] = g.Actnext[i,j,k] + born;
    return born
end

function death_event!(g::Grid, Drate, binGb, Dweight, deltat, Popvox, Necvox, K, Popgen,
    i, j, k, e)
    """
        Death event
    """
    drate = Drate*(1-binGb[1]*Dweight[1]-binGb[2]*Dweight[2]-binGb[3]*Dweight[3]);
    Pkill = (deltat/drate)*(Popvox+Necvox)/K;
    if Pkill > 1
        Pkill = 1;
    end
    if Pkill < 0
        Pkill = 0;
    end
    dead = rand(Binomial(Int64(Popgen), Pkill));
    g.Gnext[i, j, k, e] = g.Gnext[i, j, k, e] - dead;
    g.Necnext[i, j, k] = g.Necnext[i, j, k] + dead;
    return dead
end

function migration_event!(g::Grid, Migrate, binGb, Migweight, deltat, Popvox, Popgen, Necvox, K, i, j, k, e)
    """
        Migration event
    """
    migrate = Migrate*(1-binGb[1]*Migweight[1]-binGb[2]*Migweight[2]-binGb[3]*Migweight[3]);
    Pmig = (deltat/migrate)*(Popvox+Necvox)/K;
    if Pmig > 1
        Pmig = 1;
    end
    if Pmig < 0
        Pmig = 0;
    end
    migrants = rand(Binomial(Int64(Popgen),Pmig));

    neigh = 0;
    # left = migrants;
    moore = 26;
    vonN = 6;
    multinom = Multinomial(migrants,wcube);
    gone = rand(multinom);
    for movi in [-1,0,1]
        for movj in [-1,0,1]
            for movk = [-1,0,1]
                xmov = i+movi;
                ymov = j+movj;
                zmov = k+movk;
                if xmov < N+1 && ymov < N+1 && zmov < N+1 && xmov > 0 && ymov > 0 && zmov > 0 && abs(movi)+abs(movj)+abs(movk)!=0
                    neigh = neigh + 1;
                    g.Gnext[xmov, ymov, zmov, e] = g.Gnext[xmov, ymov, zmov, e] + gone[neigh];
                    g.Gnext[i, j, k, e] = g.Gnext[i, j, k, e] - gone[neigh];
                end
            end
        end
    end
end

function mutation_event!(g::Grid, Mutrate, binGb, Mutweight, deltat, Popgen, K)
    """
        Mutation event.
    """
    mutrate = Mutrate*(1-binGb[1]*Mutweight[1]-binGb[2]*Mutweight[2]-binGb[3]*Mutweight[3]);
    Pmut = (deltat/mutrate)*(Popgen/K);
    if Pmut > 1
        Pmut = 1;
    end
    if Pmut < 0
        Pmut = 0;
    end
    r = rand(1);
    r = r[1];
    if r < Pmut && e != 2^alt
        # Pick a random empty slot and turn it to mutated
        nonalter = findall(x -> x < 1, binGb);
        r2 = rand(1:length(nonalter));
        mutating = nonalter[r2];
        binGb[mutating] = 1;

        # Switch binary array back to binary string
        binGc = string(Int(binGb[1]),Int(binGb[2]),Int(binGb[3]));

        # Code below retrieves back decimal number from binary string
        decG = parse(Int, binGc, base=2)+1;
        g.Gnext[i, j, k, e] = g.Gnext[i, j, k, e] - 1;
        g.Gnext[i, j, k, decG] = g.Gnext[i, j, k, decG] + 1;
    end
end

function save_gen_space(g::Grid, t::Float64, N::Int64, subdir="files/")
    dir_to_save = joinpath(@__DIR__, subdir)
    filename = joinpath(dir_to_save, string("Gen_space_",
    string(Int64(floor(t))), ".txt"));
    open(filename, "w") do file
        for i in 1:N
            for j in 1:N
                for k in 1:N
                    if sum(g.G[i, j, k, :]) > 0
                        wpop = g.G[i, j, k, :]
                        actF = g.Act[i, j, k]
                        necF = g.Nec[i, j, k]
                        println(file, i, " ", j, " ", k, " ", wpop[1], " ",
                        wpop[2], " ", wpop[3], " ", wpop[4], " ", wpop[5],
                        " ", wpop[6], " ", wpop[7], " ", wpop[8], " ",
                        actF, " ", necF)
                    end
                end
            end
        end
    end
end

function print_curr_stats(m::Monitor, t, elapsed, evalstep)
    println("Cell no: ", m.totpop[evalstep + 1], "; Volume: ", m.Rvol[evalstep +
    1], "; Activity: ", m.Rtotnew[evalstep + 1], "; Necrotics: ",
    m.totnec[evalstep + 1], "; Het: ", m.Shannon[evalstep + 1])
    println("     Volume (alt): ", m.Vol2[evalstep + 1])
    println("Time step: ", t, "; Time elapsed: ", elapsed)
end
