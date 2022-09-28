#Problem Set 3 Code
using Plots, Parameters
@with_kw struct Primitives
    age_retire::Int64 = 46
    theta::Float64 = 0.11
    gamma::Float64 = 0.42
    sigma::Int64 = 2
    z_prod::Array{Float64,1} = [3.0, 0.5]
    birth_distribution::Array{Float64,1} = [.2037, .7963]
    markov::Array{Float64,2} = [0.9261 .0739;0.9811 0.0189] #Instantiating the asymmetric case
    alpha::Float64 = 0.36 #capital share
    delta::Float64 = 0.06 #Capital Depreciation
    beta::Float64 = 0.97
    N::Int64 = 66
    Assets::Array{Float64,1} = collect(range(0,length=250,stop=100)) ##
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




end


mutable struct Results
    pol_func:: Array{Float64,3}
    val_func:: Array{Float64,3}
    #asset_func:: Array{Float64,3}
    lab_func:: Array{Float64,3}
    v_final::Array{Float64,3}
    age_j:: Int64
    r:: Float64
    w:: Float64
    b:: Float64
    mu:: Array{Float64,2}
    K:: Float64
    L:: Float64

    #K::Array{Float64,1}
    #L::Array{Float74,1}
end

function Initialize()
    prim = Primitives()
    pol_func = zeros(prim.na,prim.nz,prim.N)
    val_func = zeros(prim.na,prim.nz,prim.N)
    lab_func = zeros(prim.na,prim.nz,prim.N)
    v_final = zeros(prim.na,prim.nz,prim.N)
    age_j = 66
    r = 0.05
    w = 1.05
    b = 0.2
    mu = zeros(prim.na,prim.nz)
    K = 3
    L = 0.5
    res = Results(pol_func,val_func, lab_func,v_final,age_j, r,w,b,mu, K, L)
    prim, res
end

function productivity(z,j)
    return z*j
end

function Backinduct(prim::Primitives,res::Results)
    @unpack age_retire, theta, gamma, sigma, z_prod, birth_distribution = prim
    @unpack markov, alpha, delta, beta, N, Assets, na, nz, age_ef = prim
    res.v_final= zeros(na,nz,N)
    w = 1.05
    r = 0.05
    b = 0.2
    K = 3.3
    L = 0.3
    v_0=1e-10
    for age_index = 66:-1:1
        v_0 = -1e10
        #if age_index == 66
        for a_index = 1:na, z_index=1:nz
            #for a_index = 1:na, z_index=1:nz
            if age_index == 66
                budget = (1+r)*Assets[a_index]
                c = budget
                if c > 0
                    val =  (c^((1-sigma)*gamma))/(1-sigma)
                    if val > v_0
                        res.pol_func[a_index,:,age_index] .= Assets[a_index]
                        res.lab_func[a_index,:,age_index] .= zero(UInt32)
                         #I am addi
                        v_0 = val
                        res.val_func[:,:,66] .= v_0 ##Add the value function over all asset levels for those in the final period of life.
                        res.v_final[:,:,66] .= v_0
                    end
                end


            elseif age_index >= 46 && age_index < 66
            #for a_index = 1:na, z_index = 1:nz
                budget = (1+r)*a + b
                for ap_index = 1:na, zp_index = 1:nz
                    c = budget - Assets[ap_index]
                    if c > 0
                        val =(c^((1-sigma)*gamma))/(1-sigma)+ beta*res.val_func(Assets[ap_index],:,age_index+1)
                        if val > v_0
                            res.pol_func[a_index,z_index,age_index] = Asset[ap_index]
                            res.lab_func[a_index,:, age_index] .= zero(UInt32)
                            v_0 = val
                        end
                    end
                end
                res.val_func[:,:,age_index] .= val
                res.v_final[:,:,age_index] .= val
        #else age_index < 46
            else age_index < 46
            #for a_index = 1:na, z_index = 1:nz
                for ap_index = 1:na, zp_index = 1:nz
                    l = (gamma*(1-theta)*productivity(z_prod[z_index],age_ef,age_index)*w-(1-gamma)*((1+r)*Assets[a_index]-Assets[ap_index]))/((1-theta)*w*productivity(z_prod[z_index],age_ef[age_index]))
                    budget = (1+r)*Assets[a_index] + w*(1-theta)*productivity(z_prod[z_index],age_ef[age_index])*l
                    c = budget - Assets[ap_index]
                    if c > 0
                        val = (((c^(gamma))*(1-l)^(1-gamma))^(1-sigma))/(1-sigma) + beta*sum(res.val_func[Assets[ap_index],:,age_index].*markov[z_index,:])
                        if val > v_0
                            res.pol_func[a_index,z_index,age_index] = Asset[ap_index]
                            res.lab_func[a_index,z_index,age_index] .= l(UInt32)
                            v_0 = val
                        end
                    end
                end
                res.val_func[:,:,age_index] .= v_0

            end
            res.v_final[:,:,age_index] .= v_0
        end
    end
    return res.v_final, res.pol_func, res.lab_func
end

prim, res = Initialize()
res.v_final= zeros(prim.na,prim.nz,prim.N)
res.v_final, res.pol_func, res.lab_func .= Backinduct(prim,res)
