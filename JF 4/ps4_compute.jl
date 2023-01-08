using DataFrames, Random, Parameters, Distributions, Accessors, CSV, LinearAlgebra
include("ps4_code.jl")

P,R = initialize(S)

overall_iteration(P,R,S)
R.exp_val_func,R.exp_val_vector = expected_val_func(P,R) 


### P_Hat Using Simulate Data
sim_phat = CCP(P,R)
R.exp_val_func_ccp,payoff_0,payoff_1,F = vbar_ccp(P,R,S,sim_phat)
#### Find Expected Value Function From CCP With Iteration Starting from Phat = 0.5

overall_ccp_iteration(P,R)
R.exp_val_func_ccp

S
