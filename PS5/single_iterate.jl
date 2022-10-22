using Parameters, Plots, Optim, Interpolations, GLM, Printf, DataFrames #last two are needed for the helpful functions
using Random, Distributions, LinearAlgebra
# cd("/Users/jacobbills/Desktop/Economics/Econ 899/PS 5/")
if occursin("Rafeh",pwd()) cd("C:/Users/Rafeh/Documents/GitHub/ECON899/PS5") end #Insert your own paths here
# cd("/Users/jacobbills/Desktop/Economics/Econ 899/PS 5/") #Insert your own paths here
if occursin("mcket",pwd()) cd("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/899/ECON899/PS5") end #Insert your own paths here

include("HelpfulFunctions_edits.jl")
include("ps5_model.jl")

P,G,S,R = initialize_overall()
Epsilon, Z = DrawShocks(S,P.N,P.T)
R.pf_k, R.pf_v = solve_HH_problem(P,G,S,R)
# R.Kappa,R.V_matrix = SimulateCapitalPath(R, P,G, Epsilon,Z)
# R.ahat0, R.ahat1, R.bhat0, R.bhat1, R.R2 = EstimateRegression(R,R.Kappa)
#First step is to draw shocks
# Solve_KS(P,G,S,R)

lambda::Float64=0.5
tol_up = 1e-3
m_gof::Float64 =1 - 1e-2
kill = 10000

@unpack a0, b0, a1, b1, R2, pf_k, pf_v = R
@unpack n_k, n_eps, n_K, n_z = G
pf_k_temp = zeros(n_k, n_eps, n_K,n_z)
pf_v_temp = zeros(n_k, n_eps, n_K,n_z)
#First step is to draw shocks
Epsilon, Z = DrawShocks(S,P.N,P.T) #use the formula from the helpful functions to get our shocks
#Second step is to set initial a0,b0,a1,b1

# R.a0, R.b0, R.a1, R.b1, R.R2, R.pf_k, R.pf_v = initialize_coef(R,P) #I have no idea how initialization works

stop = 0
count = 0
count +=1
pf_k_temp, pf_pf_v = solve_HH_problem(P,G,S,R) #note this returns K' first, different than how we normally do it
# kappa, V_b = SimulateCapitalPath(R,P,G,Epsilon,Z) #I think this is all correct now

@unpack pf_k,pf_v = R
@unpack N, T, burn = P
@unpack K_grid,k_grid = G

kappa::Array{Float64,2} = zeros(T,3) #since the first one is Kbar2. I think this is the proper bounds but I am not sure
V::Array{Float64,2} = zeros(N,T) #Storing the panel. I hate that we have to do this column by column but rows come first
pf_k_interpol = interpolate(pf_k, BSpline(Linear()) )
#do initial column seperately
for n in 1:N
    #Need to get the index for k^t_n and Kbar_t, not directly insert them.
    #V[n,1] = pf_k[K_index,trunc(Int64,Epsilon[n,1]),K_index,trunc(Int64,Z[1])] #11.55 is a hardfix, this probably doesn't work? Need to understand what the policy function is saying better
    V[n,1] = 11.556444894066576
end
kappa[1,1] = 11.556444894066576
kappa[1,2] = Z[1]
kappa[1,3] = 11.556444894066576 #Column three would be capital today and thus X?
K_val = 11.556444894066576

for t in 2:T #I think this is right because in the pseudocode, it only goes to Kbar_T
    K_index = get_index(kappa[t-1,1],K_grid)
    for n in 1:N
        #See details above about why this doesn't work right
        #V[n,t] = pf_k[trunc(Int64,get_index(V[n,t-1],K_grid)),trunc(Int64,Epsilon[n,t]),trunc(Int64,get_index(kappa[t-1,1],K_grid)),trunc(Int64,Z[t])] #need to flip my epsilon I think
        V[n,t] = pf_k_interpol(get_index(V[n,t-1],k_grid),Epsilon[n,t],K_index,Z[t-1])
    end
    kappa[t,1] = sum(V[:,t])/N
    kappa[t,2] = Z[t] #for the sort later on, tracks the state of the world when the decision was made
    kappa[t,3] = kappa[t-1,1] #This makes column three yesterday's capital choice, ie today's capital. This is X in that case (I think)
    if mod(t,1000) == 0
        println("another 10%")
    end
end

#burn the first 1000 columns/rows (kappa is weird)
kappa_b::Array{Float64,2} = kappa[burn+1:T-1,:]
V_b::Array{Float64,2} = V[:,burn+1:T-1]


R.ahat0, R.ahat1, R.bhat0, R.bhat1, R.R2 = EstimateRegression(R,R.Kappa) #Still need to write this

if (abs(R.ahat0-R.a0)+abs(R.ahat1-R.a1)+abs(R.bhat0-R.b0)+abs(R.bhat1-R.b1)) > tol_up
    R.a0, R.b0, R.a1, R.b1 = Update_Coef(R::Results, lambda) #Still need to write this up
else
    stop = 1
    println("The coefficients have converged")
    println("Goodness of fit (good, bad): ", R.R2)
end
println("a0: ",R.a0,"\tb0: ", R.b0,"\ta1: ", R.a1,"\tb1: ", R.b1, "\t R2: ", R.R2, "gof: ",abs(R.ahat0-R.a0)+abs(R.ahat1-R.a1)+abs(R.bhat0-R.b0)+abs(R.bhat1-R.b1))

println("Converged in: ",count," iterations")
println("Predicted a0: ",R.a0)
println("Predicted a1: ", R.a1)
println("Predicted b0: ", R.b0)
println("Predicted b1: ", R.b1)