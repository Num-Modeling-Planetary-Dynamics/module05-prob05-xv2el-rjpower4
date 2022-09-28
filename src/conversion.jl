# ========================================================================================
# File: conversion.jl
# Brief: Code for converting from Keplerian elements to cartesian state
# Author: Rolfe Power <rpower@purdue.edu>
# ========================================================================================

using CSV
using DataFrames

include("path_utils.jl")
include("keplerian_elements.jl")

function main()
    dataset_files = get_data_files()
    return datasets
end
