"""
    Module contains tools for testing, saving results to files, visualizastion
    etc.
"""

include("monitor.jl")

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

function reproduction_event(Popgen, Popvox, K, Necvox, Grate, deltat, i, j, k,
    e, alt, binGb)
    """
        Reproduction event.
    """
    grate = Grate*(1-binGb[1]*Gweight[1]-binGb[2]*Gweight[2]-binGb[3]*Gweight[3]);
    Prep = (deltat/grate)*(1-(Popvox+Necvox)/K);
    if Prep > 1
        Prep = 1;
    end
    if Prep < 0
        Prep = 0;
    end
    born = rand(Binomial(Int64(Popgen),Prep));
    Gnext[i,j,k,e] = Gnext[i,j,k,e] + born;
    Actnext[i,j,k] = Actnext[i,j,k] + born;
    return born
end

function death_event(Drate, binGb, Dweight, deltat, Popvox, Necvox, K, Popgen,
    i, j, k, e)
    """
        Death event.
    """
    drate = Drate*(1-binGb[1]*Dweight[1]-binGb[2]*Dweight[2]-binGb[3]*Dweight[3]);
    Pkill = (deltat/drate)*(Popvox+Necvox)/K;
    if Pkill > 1
        Pkill = 1;
    end
    if Pkill < 0
        Pkill = 0;
    end
    dead = rand(Binomial(Int64(Popgen),Pkill));
    Gnext[i,j,k,e] = Gnext[i,j,k,e] - dead;
    Necnext[i,j,k] = Necnext[i,j,k] + dead;
    return dead
end

function migration_event(Migrate, binGb, Migweight, deltat, Popvox, Popgen, Necvox, K, i, j, k, e)
    """
        Migration event.
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
                    Gnext[xmov,ymov,zmov,e] = Gnext[xmov,ymov,zmov,e] + gone[neigh];
                    Gnext[i,j,k,e] = Gnext[i,j,k,e] - gone[neigh];
                end
            end
        end
    end
end

function mutation_event(Mutrate, binGb, Mutweight, deltat, Popgen, K)
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
        Gnext[i,j,k,e] = Gnext[i,j,k,e] - 1;
        Gnext[i,j,k,decG] = Gnext[i,j,k,decG] + 1;
    end
end

function adjust_grate_migrate(Grate, Migrate, MinGrate, MaxGrate, MinMigrate,
    MaxMigrate)
    while Grate / Migrate < 0.25 || Migrate / Grate < 0.1
        Grate = rand(Uniform(MinGrate, MaxGrate));
        Migrate = rand(Uniform(MinMigrate, MaxMigrate));
    end
    Grate, Migrate
end

function save_gen_space(G::Array{Float64,4}, Act::Array{Float64,3},
    Nec::Array{Float64,3}, t::Float64, N::Int64, subdir="files/")
    dir_to_save = joinpath(@__DIR__, subdir)
    filename = joinpath(dir_to_save, string("Gen_space_",
    string(Int64(floor(t))), ".txt"));
    open(filename, "w") do file
        for i in 1:N
            for j in 1:N
                for k in 1:N
                    if sum(G[i, j, k, :]) > 0
                        wpop = G[i, j, k, :]
                        actF = Act[i, j, k]
                        necF = Nec[i, j, k]
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
