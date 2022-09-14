using Parameters, Plots, Distributed, SharedArrays #import the libraries we want
addprocs(2)
@everywhere include("/Users/jacobbills/Desktop/Economics/Econ 899/PS 1/Our Work/f22stochasticmodel_pll.jl") #import the functions that solve our growth model

@time prim, res = Initialize() #initialize primitive and results structs
@time Solve_model(prim, res) #solve the model!
@unpack val_func, pol_func = res
@unpack k_grid = prim

##############Make plots
#value function
#Plots.plot(k_grid, val_func, title="Value Function", label = ["Good Shock" "Bad Shock"])
#Plots.savefig("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/Problem Set 1/02Growthzipjulia/02_Value_Functions.png")

#policy functions
#Plots.plot(k_grid, pol_func, title="Policy Functions", label = ["Good Shock" "Bad Shock"])
#Plots.savefig("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/Problem Set 1/02Growthzipjulia/02_policy_Functions.png")

#changes in policy function
#pol_func_δ = copy(pol_func).-k_grid
#Plots.plot(k_grid, pol_func_δ, title="Policy Functions Changes", label = ["Good Shock" "Bad Shock"])
#Plots.savefig("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/Problem Set 1/02Growthzipjulia/02_polcy_changes.png")

println("All done!")
################################