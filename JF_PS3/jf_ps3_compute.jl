using Parameters, DataFrames, StatFiles, Random, Distributions, Optim, LinearAlgebra, Plots, LogExpFunctions, Distributed

#cd("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3")
#cd("/Users/jacobbills/Desktop/Economics/Econ 899/JF 3/")

ind_characteristics = DataFrame(load("./Simulated_type_distribution.dta"))
car_ch = DataFrame(load("./Car_demand_characteristics_spec1.dta"))
iv_specification = DataFrame(load("./Car_demand_iv_spec1.dta"))


#ind_characteristics = DataFrame(load("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3/PS3/Simulated_type_distribution.dta"))
#car_ch = DataFrame(load("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3/PS3/Car_demand_characteristics_spec1.dta"))
#iv_specification = DataFrame(load("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3/PS3/Car_demand_iv_spec1.dta"))

include("./ps3_code.jl")

P,R =  Initialize()
res_norm = Array{Float64,1}
res_norm = routine_overall(P,R,1,.6)
res_norm = res_norm[2:length(res_norm)]
plot(1:length(res_norm), res_norm, title="Plotting Norm Across Iterations", labels = "", legend=:topleft)
#Plots.savefig("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3/PS3/norm_plot")
Plots.savefig("./norm_plot")


##Question 2

gm_grid, min_lam, idx_lam, grid_lam = gridsearch(P,R,0,1,20,31)

plot(grid_lam, gm_grid, title="λp Grid Search", labels = "", legend=:topleft)
Plots.savefig("./gmm_obj")

println("Grid lambda: ", min_lam)


#Question  3
first_l, second_l = two_step(P,R,min_lam,31)

println("The λ found in the first step is: ", first_l)
println("The λ found in the second step is: ", second_l)