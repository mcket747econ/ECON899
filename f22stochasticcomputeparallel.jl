using Parameters, Plots #import the libraries we want
using Distributed
addprocs(2)
@everywhere include("02Growth_model.jl") #import the functions that solve our growth model

@everywehre @time prim, res = Initialize() #initialize primitive and results structs
@everywhere @time Solve_model(prim, res) #solve the model!
@everywhere @unpack val_func, pol_func = res
@everywhere @unpack k_grid = prim

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
