# ========================================================================================
# File: path_utils.jl,
# Brief: Simple functions for managing project paths
# Author: Rolfe Power <rpower@purdue.edu>
# ========================================================================================

# ----------------------------------------------------------------------------------------
# Directories Realtive to Project Root
# ----------------------------------------------------------------------------------------
project_dir() = dirname(@__DIR__)
project_dir(subs...) = joinpath(project_dir(), subs...)
data_dir(subs...) = joinpath(project_dir("data"), subs...)

# ----------------------------------------------------------------------------------------
# Getting the data files
# ----------------------------------------------------------------------------------------
function get_data_files(dir=data_dir())
    re = r"id\d{6}-XV.csv"
    output_dict = Dict{String, String}()
    candidates = readdir(dir)
    for candidate in candidates
        if match(re, candidate) !== nothing
            output_dict[candidate] = joinpath(dir, candidate)
        end
    end
    return output_dict
end
