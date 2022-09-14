#import the libraries we want

# Original Code from Fall 2021 Modified by Matt McKetty Fall 2022\
using Distributed ##Add distributed package to enable parallel processing
addprocs(2) #State the number of processes to run
@everywhere using Pkg;Pkg.add(["Plots","Parameters"]) #Add other relevant packages
@everywhere using Parameters, Plots
@everywhere include("f22stochasticmodelparallel.jl") #import the functions that solve our growth model

@everywhere @time prim, res = Initialize() #initialize primitive and results structs
@everywhere @time Solve_model(prim, res) #solve the model!
@everywhere @unpack val_func, pol_func = res #extract results for plots
@everywhere @unpack k_grid = prim #extract k array from primitive structure

##############Make plots
#value function
Plots.plot(k_grid, val_func, title="Value Function", label = ["Good Shock" "Bad Shock"])
Plots.savefig("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/Problem Set 1/02Growthzipjulia/02_Value_Functions.png")

#policy functions
Plots.plot(k_grid, pol_func, title="Policy Functions", label = ["Good Shock" "Bad Shock"])
Plots.savefig("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/Problem Set 1/02Growthzipjulia/02_policy_Functions.png")

#changes in policy function
pol_func_δ = copy(pol_func).-k_grid
Plots.plot(k_grid, pol_func_δ, title="Policy Functions Changes", label = ["Good Shock" "Bad Shock"])
Plots.savefig("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/Problem Set 1/02Growthzipjulia/02_polcy_changes.png")

println("All done!")
################################
