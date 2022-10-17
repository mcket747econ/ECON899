using Parameters, Plots, Optim, Interpolations #last two are needed for the helpful functions 
using Random
# cd("/Users/jacobbills/Desktop/Economics/Econ 899/PS 5/")
if occursin("Rafeh",pwd()) cd("C:/Users/Rafeh/Documents/GitHub/ECON899/PS5") end #Insert your own paths here 
#cd() #Insert your own paths here

include("HelpfulFunctions_edits.jl")
include("PS5_model.jl")

Initialize() #This isn't defined in either file yet, still needs to be done 
P = Params()
S = Shocks()
Solve_KS(P,G,S,P,0.5, 1e-3, 1 - 1e-2, 100) #Have it set to kill pretty fast, should be changed in the final run I guess 



