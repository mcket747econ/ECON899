Original_agrid=1; 
using DelimitedFiles
using Interpolations, Plots, Parameters, DataFrames, Random, Distributions, GLM, Optim, Printf


@with_kw struct Params
    beta::Float64=0.8;
    A::Float64=0.6;

    θ::Float64=0.64; 
    cf::Float64=12; 
    ce::Float64=40; 
    
    tau::Float64=0; 
    p::Float64=1; 
    w::Float64=1; 
    
    rho::Float64=0.93; 
    sigma_logz::Float64=sqrt(0.53); 
    sigma_epsilon::Float64=sqrt((1-rho)*((sigma_logz)^2));
    a::Float64=0.078; 

    Fs = reshape(readdlm("C:\\Users\\Rafeh\\Documents\\GitHub\\ECON899\\PS6\\Fs.txt", ',', Float64),5,5)'
    vs = readdlm("C:\\Users\\Rafeh\\Documents\\GitHub\\ECON899\\PS6\\vs.txt", ',', Float64)
    s_grid = [3.98e−4, 3.58, 6.82, 12.18, 18.79]'
    tol = 1e-3
    n_s = length(s_grid)
end

function n_opt(p::Float64,s::Float64,Params::Params)
    @unpack θ = Params
    (1/θ*p*s)^(1/θ-1) 
end 

function π(p::Float64,s::Float64,Params::Params)
    @unpack θ, cf = Params
    n_star = n_opt(p,s,Params)
    p*s*n_star - n_star - p*cf
end    

# s::Float64
function W(p::Float64,Params::Params)        ## Don't think this is correct - this would be just two periodss
    @unpack θ, n_s = Params
    W::Array{Float64,1} = zeros(0,4)
    for is = 1:n_s
        W_int = 0.0 
        for isp in 1:n_s
            W_int += Fs[is,isp]*π(p,s,Params)           ## Placeholder - need to program the DP properly
        end
        W[is] = W_int
    end
    W       ##Returns a 1x4 array
end

function W_ent(s::Float64,p::Float64,Params::Params)        ## Don't think this is correct - this would be just two periodss
    @unpack θ, n_s = Params
    W::Array{Float64,1} = 0.0
    for is = 1:n_s
        W_int += vs[is]*π(p,s,Params)
    end
    W
end

function VF(s::Float64,p::Float64,Params::Params) 
    @unpack θ, n_s = Params
    W_int::Array{Float64,1} = W(p,params)
    x = 0
    W = π(p,s)
    if W_int[s]>0
        x = 1
        W+= W+W_int 
    end
    x, W
end


