using Parameters, Plots, Optim, Interpolations, GLM, Printf, DataFrames #last two are needed for the helpful functions
using Random, Distributions, LinearAlgebra

# cd("/Users/jacobbills/Desktop/Economics/Econ 899/PS 5/")
if occursin("Rafeh",pwd()) cd("C:/Users/Rafeh/Documents/GitHub/ECON899/PS5") end #Insert your own paths here
# cd("/Users/jacobbills/Desktop/Economics/Econ 899/PS 5/") #Insert your own paths here
if occursin("mcket",pwd()) cd("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/899/ECON899/PS5") end #Insert your own paths here

include("HelpfulFunctions_edits.jl")
include("ps5_model.jl")

## Moved everything to solver
## Should probably remove "print" in KS_solve
overall_solve()


