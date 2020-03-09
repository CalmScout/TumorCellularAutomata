"""
    Module contains tools for testing, saving results to files, visualizastion
    etc.
"""

@everywhere include("monitor.jl")
@everywhere include("grid.jl")

function compare_files(path1::String, path2::String)
    """
    Reads two text files from provided pathes and compares them.
    """
    return read(path1, String) == read(path2, String)
end

function compare_folders(path_folder_1::String, path_folder_2::String)
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

function decimal2binstr(e::Int64, alt::Int64)
    """
        Converts decimal number to binary string.
    """
    binG = digits(e-1, base=2, pad=alt) |> reverse
    return convert(Array{Float64, 1}, binG)
end

function save_gen_space(g::Grid, t::Int64, N::Int64, subdir::String="files/")
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

function print_curr_stats(m::Monitor, t::Int64)
    println("Cell no: ", m.totpop[m.evalstep], "; Volume: ", m.Rvol[m.evalstep],
     "; Activity: ", m.Rtotnew[m.evalstep], "; Necrotics: ",
    m.totnec[m.evalstep], "; Het: ", m.Shannon[m.evalstep])
    println("     Volume (alt): ", m.Vol2[m.evalstep])
    println("Time step: ", t, "; Time elapsed: ", m.elapsed)
end

function get_non_zeroes(arr::SharedArray{CartesianIndex{3}, 1})
    """
    Return slice of array till first zero CartesianIndex.
    """
    function get_first_zero_idx(arr::SharedArray{CartesianIndex{3}, 1})
        result = -1
        for idx in 1:length(arr)
            if arr[idx] == CartesianIndex(0, 0, 0)
                result = idx
                break
            end
        end
        return result
    end
    return arr[1 : get_first_zero_idx(arr) - 1]
end
