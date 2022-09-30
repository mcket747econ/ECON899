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
    Assets::Array{Float64,1} = collect(range(0,length=2000,stop=50)) ##updated to be bigger
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
    productivity_e::Matrix{Float64} = age_ef*z_prod'


##Make a matrix


end


mutable struct Results
    pol_func:: Array{Float64,3}
    val_func:: Array{Float64,3}
    #asset_func:: Array{Float64,3}
    lab_func:: Array{Float64,3}
    #v_final::Array{Float64,3}
    age_j:: Int64
    r:: Float64
    w:: Float64
    b:: Float64
    mu:: Array{Float64,3} #mu should have 3 slots, (j,a,z) like everything else 
    K:: Float64
    L:: Float64
    K1:: Float64 #result of market clearing
    L1:: Float64 #result of market clearing 
    #K::Array{Float64,1}
    #L::Array{Float74,1}
    K_new:: Array{Float64,3} #defining these outside my equation
    L_new:: Array{Float64,3}
end

function Initialize()
    prim = Primitives()
    pol_func = zeros(prim.na,prim.nz,prim.N)
    val_func = zeros(prim.na,prim.nz,prim.N)
    lab_func = zeros(prim.na,prim.nz,prim.N)
    #v_final = zeros(prim.na,prim.nz,prim.N)
    age_j = 66
    r = 0.05
    w = 1.05
    b = 0.2
    mu = zeros(prim.na,prim.nz,prim.N) #updated with the third slot 
    K = 3
    L = 0.5
    K1 = 0 #doesn't really matter, get filled in later 
    L1 = 0  #doesn't really matter, gets filled in later
    K_new = zeros(prim.na,prim.nz,prim.N)
    L_new = zeros(prim.na,prim.nz,prim.N-prim.age_retire)
    res = Results(pol_func,val_func, lab_func,age_j, r,w,b,mu, K, L, K1, L1, K_new, L_new)
    prim, res
end

# function productivity(z,j)
#     return z*j
# end

function Backinduct(prim::Primitives,res::Results)
    @unpack age_retire, theta, gamma, sigma, z_prod, birth_distribution, delta = prim
    @unpack markov, alpha, delta, beta, N,Assets, na, nz, age_ef, productivity_e = prim
    @unpack w,r,b= res
    #res.v_final= zeros(na,nz,N)
    # w = 1.05
    # r = 0.05
    # b = 0.2
    # K = 3.3
    # L = 0.3

    for age_index = 66:-1:1
        #if age_index == 66

        for a_index = 1:na, z_index=1:nz
            v_0 = -1e10
            #for a_index = 1:na, z_index=1:nz
            if age_index == 66

                budget = (1+r-delta)*Assets[a_index] + b
                c = budget
                if c > 0
                    val =  (c^((1-sigma)*gamma))/(1-sigma)
                    if val > v_0
                        res.pol_func[a_index,:,age_index] .= Assets[a_index]
                        res.lab_func[a_index,:,age_index] .= zero(UInt32)
                         #I am addi
                        v_0 = val
                        res.val_func[a_index,:,age_index] .= v_0 ##Add the value function over all asset levels for those in the final period of life.
                        #res.v_final[:,:,66] .= v_0
                    end
                end


            elseif age_index >= 46 && age_index < 66

            #for a_index = 1:na, z_index = 1:nz
                budget = (1+r-delta)*Assets[a_index] + b
                for ap_index = 1:na, zp_index = 1:nz
                    c = budget - Assets[ap_index]
                    if c > 0
                        val =(c^((1-sigma)*gamma))/(1-sigma)+ beta*res.val_func[ap_index,zp_index,age_index+1]
                        if val > v_0
                            res.pol_func[a_index,:,age_index] .= Assets[ap_index]
                            res.lab_func[a_index,:, age_index] .= zero(UInt32)
                            v_0 = val
                            res.val_func[a_index,:,age_index] .= v_0
                        end
                    end
                end
                #res.val_func[:,:,age_index] .= val
                #res.v_final[:,:,age_index] .= val
        #else age_index < 46
            else age_index < 46
                if z_index == 1
                    v0_1 = -1e15
            #for a_index = 1:na, z_index = 1:nz
                    for ap_index = 1:na, #zp_index = 1:nz
                        l = 1.0*(gamma*(1-theta)*productivity_e[age_index,z_index]*w-(1-gamma)*((1+r)*Assets[a_index]-Assets[ap_index]))/((1-theta)*w*productivity_e[age_index,z_index])
                        if l < 0
                            l = 0
                        elseif l > 1
                            l = 1
                        end
                        budget = (1+r-delta)*Assets[a_index] + w*(1-theta)*productivity_e[age_index,z_index]*l
                        c = budget - Assets[ap_index]
                        if c > 0
                            val = (((c^(gamma))*(1-l)^(1-gamma))^(1-sigma))/(1-sigma) + beta*sum(res.val_func[ap_index,:,age_index+1].*markov[z_index,:])
                            if val > v0_1
                                res.pol_func[a_index,z_index,age_index] = Assets[ap_index]
                                res.lab_func[a_index,z_index,age_index] = 1.0* l
                                v0_1 = val
                                res.val_func[a_index,z_index,age_index] = 1.0* v0_1
                            end
                        end
                    end
                else
                    v0_2 = -1e20
            #for a_index = 1:na, z_index = 1:nz
                    for ap_index = 1:na, #zp_index = 1:nz
                        l = 1.0* (gamma*(1-theta)*productivity_e[age_index,z_index]*w-(1-gamma)*((1+r)*Assets[a_index]-Assets[ap_index]))/((1-theta)*w*productivity_e[age_index,z_index])
                        if l < 0
                            l = 0
                        elseif l > 1
                            l = 1
                        end
                        budget = (1+r-delta)*Assets[a_index] + w*(1-theta)*productivity_e[age_index,z_index]*1.0*l
                        c = budget - Assets[ap_index]
                        if c > 0
                            val = (((c^(gamma))*(1-l)^(1-gamma))^(1-sigma))/(1-sigma) + beta*sum(res.val_func[ap_index,:,age_index+1].*markov[z_index,:])
                            if val > v0_2
                                res.pol_func[a_index,z_index,age_index] = Assets[ap_index]
                                res.lab_func[a_index,z_index,age_index] = 1.0*l
                                v0_2 = val
                                res.val_func[a_index,z_index,age_index] = 1.0*v0_2


                            end
                        end
                    end
                #res.val_func[:,:,age_index] .= v_0 You've already taken care of this within the previous loops
                end
            end
            #res.v_final[:,:,age_index] .= v_0
        end
    end
    return res.val_func, res.pol_func, res.lab_func
end

prim, res = Initialize()
#res.v_final= zeros(prim.na,prim.nz,prim.N)
@time Backinduct(prim,res)


 ####################################Question 2



 ###############Question 3
function MarketClearing(prim::Primitives, res::Results)
    for j in 1:N #Probably a way to do this better, but I'll need to see the outputs first 
        for m in 1:na
            for z in 1:2
                 K_new[m,z,j] = F[m,z,j]*Assets[m] #F is a placeholder here for whatever we get out of the stat dist 
                if j < age_retired
                    L_new[m,z,j] = F[m,z,j]*productivity_e[j,z]*lab_func[m,z,j] #I made some assumptions on which things I was supposed to pulle  
                else
                    continue #does Julia need this? 
                end
            end
        end
    end
        K1 = sum(K_new)
        L1 = sum(L_new)
        r = (alpha*L1^(1-alpha))/K1^(1-alpha)
        w = ((1-alpha)*K1^(alpha))/L1^(alpha)
        retired_mass = 0
        for j in age_retired:N #or can do a sum across those numbers in mu, I'm not sure how this will be stored from stat distribution function
            retired_mass += (1/1.011)^j 
        end
        #retired_mass = sum(mu(age_retired:N,:,:)) #sum alternative once I know how mass is stored    
        b = (theta*w*L1)/retired_mass
    return K1, L1, r, w, b 
end

function KL_update(res::Results, lam::Float64 = .01) #Marking part of the loop a function 
    @unpack K1, L1, K, L = res        
    K = lam*K1+(1-lam)*K 
    L = lam*L1+(1-lam)*L
    return K, L
end 

function Run_all(prim::Primitives,res::Results, tol::Float64 = 1e-3)
    stop = 0
    while stop == 0
        Backinduct(prim,res)
        CalculateStatDist(prim,res) ##rename with the appropriate function
        K1, L1, r, w, b = MarketClearing(prim,res)
        if (abs(K1-K) + abs(L1-L)) > tol
            KL_update(res,.3) #can adjust as needed, .01 is what the prompt says but supposedly we can go fast
        else
            stop = 1
        end 
    end
end

