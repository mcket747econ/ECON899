#Problem Set 3 Code
using Plots, Parameters
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
    Assets::Array{Float64,1} = collect(range(0,length=500,stop=18.75)) ##made it shorter for now while bug testing, each step in .0375
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
    F:: Array{Float64,3} # Distribution over asset holdings, age, and productivity
    K:: Float64
    L:: Float64
    K1:: Float64 #result of market clearing
    L1:: Float64 #result of market clearing
    #K::Array{Float64,1}
    #L::Array{Float74,1}
    K_new:: Array{Float64,3} #defining these outside my equation
    L_new:: Array{Float64,3}
    welf:: Float64
    CV:: Float64
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
    F = zeros(prim.na,prim.nz,prim.N)
    K = 3.3
    L = 0.3
    K1 = 0 #doesn't really matter, get filled in later
    L1 = 0  #doesn't really matter, gets filled in later
    K_new = zeros(prim.na,prim.nz,prim.N)
    L_new = zeros(prim.na,prim.nz,prim.age_retire-1)
    welf = 0.0
    CV = 0.0
    res = Results(pol_func,val_func, lab_func,age_j, r,w,b,F, K, L, K1, L1, K_new, L_new, welf, CV)
    prim, res
end

function labor(prim::Primitives,res::Results,age_index,a_index,ap_index,z_index)  ##Optimal Labor Supply function
    @unpack theta, gamma, produc,Assets = prim
    @unpack w,r = res
    l = (gamma*(1-theta)*produc[age_index,z_index]*w-(1-gamma)*((1+r)*Assets[a_index]-Assets[ap_index]))/((1-theta)*w*produc[age_index, z_index]) ##Optimal Endogenous labor supply, including depreciation
    if l < 0
        l = 0
    elseif l > 1
        l = 1
    end
    l
end

function utility(prim::Primitives, res::Results, c, l,age_index)  ##Utility functions
    @unpack sigma, gamma = prim
    u = -Inf
    if  age_index > 45  ##Utility functions for retired population
        if c>0
            u = (c^((1-sigma) * gamma))/(1 - sigma)
        end
    elseif age_index < 46 #Utility function for working population
        if c > 0
            u = (((c^gamma) * ((1 - l)^(1-gamma)))^(1-sigma))/(1 - sigma)
        end
    end

    u
end

function Backinduct(prim::Primitives,res::Results)  #Function to iterate backwards over ages to determine labor supply, asset holdings, as a function of producitivity, and age and future period asset holding
    @unpack age_retire, theta, gamma, sigma, z_prod, birth_distribution = prim
    @unpack markov, alpha, delta, beta, N, Assets, na, nz, age_ef, produc = prim
    @unpack w, r, b = res
    # res.v_final= zeros(na,nz,N)
    # w = 1.05 #######
    # r = 0.05 #######
    # b = 0.2 #######
    # K = 3.3 #######
    # L = 0.3 #######
    # v_0=1e-10
    for age_index = 66:-1:1
        #println("Age: ", age_index)
        #if age_index == 66
        for a_index = 1:na, z_index=1:nz #Iterate over potential asset holdings

            #for a_index = 1:na, z_index=1:nz
            if age_index == 66   ##The problem for those in the last period of life
                res.val_func[a_index,z_index, age_index] = -Inf #Set a very low initial guess for value function
                budget = (1+r)*Assets[a_index] + b #
                c = budget
                l = 0
                val = utility(prim, res,c,l,age_index) #Our value function in the last period of life involves no saving
                #val =  (c^((1-sigma)*gamma))/(1-sigma)
                if val > res.val_func[a_index,z_index, age_index]   ##Check if our value function is greater than initital guess
                    res.val_func[a_index,z_index,age_index] = val #update value of this value function
                    res.pol_func[a_index,:,age_index] .= zero(UInt32) #We set both labor and saving to zero in this period
                    res.lab_func[a_index,:,age_index] .= zero(UInt32)
                     #I am addi
                    ##Add the value function over all asset levels for those in the final period of life.
                    # res.v_final[:,:,66] .= v_0
                end



            elseif age_index >= 46 && age_index < 66  ##The problem for retired agents
                res.val_func[a_index,z_index, age_index] = -Inf #Initial value function
                l = 0
            #for a_index = 1:na, z_index = 1:nz
                budget = (1+r)*Assets[a_index] + b
                for ap_index = 1:na, zp_index = 1:nz  ##We now care about saving. Include state just for consistency
                    c = budget - Assets[ap_index]

                    val =utility(prim, res,c,l,age_index)+ beta*res.val_func[ap_index,zp_index,age_index+1]
                    if val > res.val_func[a_index,z_index,age_index]   ##Conduct the same progressions as above. Updating the value function if a level of saving produces higher continuation value than the previous
                        res.val_func[a_index,z_index, age_index] = val ##level of saving
                        res.pol_func[a_index,:, age_index] .= Assets[ap_index]  ##we update our saving function
                        res.lab_func[a_index,:, age_index] .= zero(UInt32) #We dont work
                         ##Add the value function over all asset levels for those in the final period of life.
                    end

                end
                # res.val_func[:,:,age_index] .= val
                # res.v_final[:,:,age_index] .= val
        #else age_index < 46
    else age_index < 46  ##the problem for workers
                res.val_func[a_index,z_index, age_index] = -Inf
            #for a_index = 1:na, z_index = 1:nz
                for ap_index = 1:na
                    l = labor(prim,res,age_index,a_index,ap_index,z_index)

                    budget = (1+r)*Assets[a_index] + w*(1-theta)*produc[age_index,z_index]*l  ##Our state is now very important
                    c = budget - Assets[ap_index]
                    val = utility(prim, res,c,l,age_index) + beta*sum(res.val_func[ap_index,:,age_index+1].*markov[z_index,:]) ##value function now heavily impacted by expectation over future productivity states
                    if val > res.val_func[a_index,z_index,age_index] #Same progression as before
                        res.val_func[a_index,z_index, age_index] = val
                        res.pol_func[a_index,z_index,age_index] = Assets[ap_index] #Update our saving
                        res.lab_func[a_index,z_index,age_index] = l ##Now we have to set our
                        ##Add the value function over all asset levels for those in the final period of life.
                    end
                end
                # res.val_func[:,:,age_index] .= v_0

            end
            # res.v_final[:,:,age_index] .= v_0
        end
    end
    return res.val_func, res.pol_func, res.lab_func
end


#prim, res = Initialize()
#res.v_final= zeros(prim.na,prim.nz,prim.N)
#@time Backinduct(prim,res)

##Note, make sure to run Backinduct with the non-endogenous L,K,r,w,b for the plots and question 1


 ####################################Question 2
 function Fdist_sum(prim::Primitives,res::Results)  ##Function to calculate distribution
    @unpack N,n,nz,na, Assets, markov=prim
    @unpack F = res
    F = zeros(na,nz,N) #reinitialization each loops fixes weird distribution being too big problem
    pop_weights = ones(66) #Initial weights
    #mu = zeros(N,nz)
    for j=1:N-1
        pop_weights[j+1] = pop_weights[j]/(n+1)
    end  ##Instantiate relative population sizes
    pop_sum = sum(pop_weights)
    pop_normal = pop_weights/pop_sum #Normalize the population
    F[1,1,1] = .2037 #weights of population with initital productivities
    F[1,2,1] = .7963
    # s=-Inf
    # s_0 = -Inf
    # for j=1:66
    #     s += (1/(1+.011))^(j)
    #     if s >= s_0
    #         s_0=s
    #     end
    #     s_0
    # end
    # s_0
    for age_index = 1:65  #for all ages
      for a_index = 1:na, z_index = 1:nz #Over all asset levels and productivity levels
          if F[a_index, z_index, age_index] > 0 #If there is a nonzero measure of the population with a given asset and productivity level
              ap = res.pol_func[a_index, z_index,age_index] #Find the policy function prescription for savings given a persons age, holdings and state
              ap_index = argmin(abs.(ap.-Assets)) #find the assets prescribed by policy function on the asset grid
              for zp_index= 1:nz
                  F[ap_index,zp_index,age_index + 1 ] += markov[z_index, zp_index]*F[a_index, z_index, age_index] ##Given the measure of the population with a given level of assets, productivity and age
                  #find the measure of the population that will end up with the prescribed policy function level of assets in the next year of life by weighting by the markov matrix
              end
         end
     end
 end
 for a_index= 1:na,z_index = 1:nz  ###Multiply each measure of the distribution by the normalized population weights
     F[a_index,z_index,:] = F[a_index,z_index,:].*pop_normal
 end

    return pop_weights,pop_sum,pop_normal, F
end



 ###############Question 3
function MarketClearing(prim::Primitives, res::Results) #check these, might be overcounting something?
    @unpack N, na, nz, alpha, age_retire, Assets, produc, theta = prim
    @unpack K1, L1, L, K, r, w, b, lab_func, F, K_new, L_new = res
    for j in 1:N #Probably a way to do this better, but I'll need to see the outputs first, should be N but it didn't like that for some reason
        for m in 1:na
            for z in 1:2
                 res.K_new[m,z,j] = F[m,z,j]*Assets[m] #F is a placeholder here for whatever we get out of the stat dist
                if j < age_retire
                    res.L_new[m,z,j] = F[m,z,j]*produc[j,z]*lab_func[m,z,j] #I made some assumptions on which things I was supposed to pulle
                else
                    continue #does Julia need this?
                end
            end
        end
    end
        res.K1 = sum(K_new)
        res.L1 = sum(L_new)
    return res.K1, res.L1
end

function KL_update(prim::Primitives,res::Results, pop_normal, lam::Float64 = .01) #Marking part of the loop a function
    @unpack alpha, theta, age_retire, N, delta = prim
    @unpack K1, L1, K, L, w, r, b = res
    res.K = lam*K1+(1-lam) *K
    res.L = lam*L1+(1-lam)*L
    res.r = (alpha*(L^(1-alpha)))/(K^(1-alpha)) - delta
    res.w = ((1-alpha)*(K^(alpha)))/(L^(alpha))
    retired_mass = sum(pop_normal[age_retire:N])
    #for j in age_retired:N #or can do a sum across those numbers in mu, I'm not sure how this will be stored from stat distribution function
     #   retired_mass += (1/1.011)^j
    #end
    #retired_mass = sum(mu(age_retired:N,:,:)) #sum alternative once I know how mass is stored
    res.b = (theta*w*L)/retired_mass
    return res.K, res.L, res.r, res.w, res.b
end

function welfare_fctn(prim::Primitives, res::Results)
    ##Get relevant parameters for calculation; asset grid, length of asset grid,
    ##number of states and number of time periods
    @unpack Assets, na, nz, N = prim
    a_grid = zeros(na,nz, N)       ##Initialize grid to length Time x Assets x States
    welfare = res.val_func .* res.F     ##welfare grid - Value Function x distribution
    welfare = welfare[isfinite.(welfare)]
    res.welf = sum(welfare)

    for j=1:N, z_index=1:nz
        a_grid[ :, z_index,j] = Assets
    end

    #mean_welfare = Base.mean(a_grid,res.F)
    mean_welfare =sum(a_grid.*res.F)      ##Weight the assets by the distribution
    var_welfare = sum(res.F.* (a_grid.^2)) - mean_welfare.^2
    res.CV = sqrt(var_welfare)/mean_welfare  
    res
end

#welfare_fctn(prim,res)


################
function Run_all(prim::Primitives,res::Results, tol::Float64 = 1e-3, lam::Float64 = .01, kill::Int64 = 1000)
    res.r = (prim.alpha*(res.L^(1-prim.alpha)))/(res.K^(1-prim.alpha)) - prim.delta
    res.w = ((1-prim.alpha)*(res.K^(prim.alpha)))/(res.L^(prim.alpha))
    retired_mass = 0
    for j in prim.age_retire:prim.N
        retired_mass += (1/1.011)^j
    end
    res.b = (prim.theta*res.w*res.L)/retired_mass
    stop = 0
    count = 1
    while  count <= kill && stop == 0
        #print("iteration: ", count, "     ")
        Backinduct(prim,res)
        pop_Weights, pop_sum, pop_normal, res.F = Fdist_sum(prim,res) ##rename with the appropriate function
        K1, L1 = MarketClearing(prim,res)
        if (abs(K1-res.K) + abs(L1-res.L)) > tol
            #print((abs(K1-res.K) + abs(L1-res.L)), "       ")
            K, L, r, w, b = KL_update(prim, res, pop_normal, lam) #can adjust as needed, .01 is what the prompt says but supposedly we can go fast
            count += 1
        else
            stop = 1
        end
    end
    welfare_fctn(prim,res) #get aggregate welfare
    print("Finished in ",count," iterations")
    println("Capital: ", res.K)
    println("Labor: ", res.L)
    println("Rent: ", res.r)
    println("wages: ", res.w)
    println("b: ", res.b)
    println("Welfare: ", res.welf)
    println("CV: ", res. CV)
end

###The real test
prim, res = Initialize()
pop_weights,pop_sum,pop_normal, F =  Fdist_sum(prim,res)
@time Run_all(prim,res, 1e-1, .12)
