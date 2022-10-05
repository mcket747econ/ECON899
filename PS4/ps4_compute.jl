using Parameters, Plots, Printf, Setfield, DataFrames, DelimitedFiles
cd("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/899/ECON899/PS4")
# Pkg.add("Accessors")
include("ps_4_code.jl")
prim, res = Initialize()

val_func_0, K_0, L_0, r_0, w_0, b_0, F_0 =  Run_all(prim,res, .01, .3)
prim_2 = @set prim.theta = 0
val_func_n, K_n, L_n, r_n, w_n, b_n, F_n = Run_all(prim_2,res, .01, .3)


mutable struct new_primitives
    T_delta::Float64
    K_t::Array{Float64,1}
    L_t::Array{Float64,1}








end

function Initialize2()
    T_delta = (K_n - K_0 )/prim.T
    K_t[t] = K_0 + prim.t_Array[t]*delta
    L_t::Array{Float64,1}
    res2 = new_primitives(T_delta, K_t, L_t)







end
res2 = Initialize2

function shoot_backward()
    @unpack T_delta, K_t, L_t, res2 = res2
    @unpack na = prim
    for t in T-1:-1:1
        for a = 1:na
            for z = 1:nz
                res.K_t[t] = K_0 + t*T_delta
                res.val_func[a,z,t] =  







end








end





###Now We Compute the Transition Path


prim, res = Initialize()
elapse =






@with_kw struct Primitives
    age_retire::Int64 = 46
    theta::Float64 = 0.11
    gamma::Float64 = .42
    sigma::Int64 = 2
    z_prod::Array{Float64,1} = [3.0, 0.5]
    birth_distribution::Array{Float64,1} = [.2037, .7963]
    markov::Array{Float64,2} = [0.9261 .0739;0.0189 0.9811 ] #Instantiating the asymmetric case
    alpha::Float64 = 0.36 #capital share
    delta::Float64 = 0.06 #Capital Depreciation
    beta::Float64 = 0.97
    N::Int64 = 66
    Assets::Array{Float64,1} = collect(range(0,length=1000,stop=30)) ##made it shorter for now while bug testing, each step in .0375
    na::Int64 = length(Assets)
    nz::Int64 = length(z_prod)
    age_ef::Array{Float64,1} = [0.59923239,0.63885106, .67846973,
    0.71808840, .75699959, 0.79591079, 0.83482198,
      0.87373318, 0.91264437,0.95155556, 0.99046676,
       0.99872065, 1.0069745, 1.0152284, 1.0234823,
       1.0317362, 1.0399901, 1.0482440, 1.0564979,
       1.0647518, 1.0730057, 1.0787834, 1.0845611,
       1.0903388, 1.0961165, 1.1018943, 1.1076720,
       1.1134497, 1.1192274, 1.1250052, 1.1307829,
       1.1233544, 1.1159259, 1.1084974, 1.1010689,
       1.0936404,
       1.0862119,
       1.0787834,
       1.0713549,
       1.0639264,
       1.0519200,
       1.0430000,
       1.0363000,
       1.0200000,
       1.0110000]
    produc = age_ef * z_prod' ##Productivity at different life stages
    n::Float64 = .011
##Make a matrix

end
