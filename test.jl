"""
    Module for testing modifications in the code during refactoring.
"""

include("tools.jl")

# compare folders `files` after refactoing and folder `files_test`
path_orig = joinpath(@__DIR__, "files_test/")
path_refactored = joinpath(@__DIR__, "files/")
println("Result of folders content testing: ",
compare_folders(path_orig, path_refactored))
