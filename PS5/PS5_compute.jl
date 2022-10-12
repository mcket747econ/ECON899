using Parameters, Plots, Optim, Interpolations #last two are needed for the helpful functions 

cd("/Users/jacobbills/Desktop/Economics/Econ 899/PS 5/")
#cd() #Insert your own paths here 
#cd() #Insert your own paths here

include("PS5_model.jl")
include("HelpfulFunctions_edits.jl")

Initialize() #This isn't defined in either file yet, still needs to be done 

Solve_KS(P,G,S,P,0.5, 1e-3, 1 - 1e-2, 100) #Have it set to kill pretty fast, should be changed in the final run I guess 



