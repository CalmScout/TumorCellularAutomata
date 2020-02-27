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

    binG = string(e-1, base=2)     # First of all, we retrieve binary string
    while length.(binG) < 3         # As the string may not be of length 0, we solve this
        binG = string("0", binG)
    end
    binGa = split(binG, "")         # Now we need an array instead of a string
    binGb = zeros(length(binGa))   # We create a new variable to store int array
    for i in 1:length(binGa)        # We convert array elements from char to int
        binGb[i] = parse(Int, binGa[i])
    end
    return binGb
end

function save_gen_space(g::Grid, t::Float64, N::Int64, subdir="files/")
    dir_to_save = joinpath(@__DIR__, subdir)
    filename = joinpath(dir_to_save, string("Gen_space_",
    string(Int64(floor(t))), ".txt"))
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

function print_curr_stats(m::Monitor, t)
    println("Cell no: ", m.totpop[m.evalstep], "; Volume: ", m.Rvol[m.evalstep],
     "; Activity: ", m.Rtotnew[m.evalstep], "; Necrotics: ",
    m.totnec[m.evalstep], "; Het: ", m.Shannon[m.evalstep])
    println("     Volume (alt): ", m.Vol2[m.evalstep])
    println("Time step: ", t, "; Time elapsed: ", m.elapsed)
end
