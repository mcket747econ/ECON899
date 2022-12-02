using Parameters, DataFrames, StatFiles, Random, Distributions, Optim, LinearAlgebra, Plots
cd("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3")

ind_characteristics = DataFrame(load("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3/PS3/Simulated_type_distribution.dta"))
car_ch = DataFrame(load("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3/PS3/Car_demand_characteristics_spec1.dta"))
iv_specification = DataFrame(load("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3/PS3/Car_demand_iv_spec1.dta"))

include("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3/jf_ps3_code.jl")

P,R =  Initialize()
res_norm = Array{Float64,1}
res_norm = routine_overall(P,R,1)
res_norm = res_norm[2:length(res_norm)]
plot(1:length(res_norm), res_norm, title="Plotting Norm Across Iterations", labels = "", legend=:topleft)
Plots.savefig("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3/PS3/norm_plot")
