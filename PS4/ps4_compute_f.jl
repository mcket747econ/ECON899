using Parameters, Plots, Printf, Setfield, DataFrames, DelimitedFiles
cd("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/899/ECON899/PS4")
# Pkg.add("Accessors")
include("ps_4_code_f.jl")
prim, res = Initialize()

val_func_0,pol_func_0,lab_func_0, K_0, L_0, r_0, w_0, b_0, F_0 =  Run_all(prim,res, .01, .3) ##Run code from pset3 to find initital steady state values with social security
prim_2 = @set prim.theta = 0
val_func_n, pol_func_n,lab_func_n, K_n, L_n, r_n, w_n, b_n, F_n = Run_all(prim_2,res, .01, .3) ##Running code from pset3 to find steady state without social security


mutable struct new_primitives  ##Initializing new mutable structure
    #T_delta::Float64
    #l_delta::Float64
    K_t::Array{Float64,1}      ##Transition values of aggregate capital
    K_potential::Array{Float64,1} #Initializing array that will be used later
    L_t::Array{Float64,1} ##aggregate labor transition array
    # val_func_n::Array{Float64,1}
    sav_func_t::Array{Float64,4}
    val_func_t::Array{Float64,4}#value function transition
    lab_func_t::Array{Float64,4} #labor function transition
    pol_func_t::Array{Float64,4} #policy function transition
    F_t::Array{Float64,4}
    #T::Int64
    #T_array::Array{Int64,1}
    pop_normal::Array{Float64,1}
    Assets::Array{Float64,4}
    w_t::Array{Float64, 1}                                                    # transition path of wage
    r_t::Array{Float64, 1}                                                    # transition path of interest rate
    b_t::Array{Float64, 1}
    theta::Array{Float64} #array of capital tax values
end


function Initialize2(T,t_iter)   ###Initialize our main structures
    #T = 30
    #t_iter = 1
    #T_array = collect(range(1,length=30,stop=T))
    #T_delta = (K_n - K_0 )/prim.T
    #l_delta = (L_n - L_0)/prim.T
    K_t = collect(range(K_0, stop =K_n, length = T))  #initial linear guess of capital transition path
    K_potential = zeros(T)
    L_t = collect(range(L_0, stop = L_n, length = T))
    # K_t_init = collect(range(K_0, stop =K_n, length = T))
    # L_t_init = collect(range(L_0, stop = L_n, length = T))
    val_func_t = zeros(prim.na,prim.nz,prim.N,T) ###initialize various functions as
    sav_func_t = zeros(prim.na,prim.nz,prim.N,T)
    lab_func_t = zeros(prim.na,prim.nz,prim.N,T)
    pol_func_t = zeros(prim.na,prim.nz,prim.N,T)
    F_t = zeros(prim.na,prim.nz,prim.N,T) #Initialize distribution function
    F_t[:,:,:,1] = F_0
    pop_normal = ones(66)
    Assets = zeros(prim.na,prim.nz,prim.N,T)
    r_t = (prim.alpha.*(L_t.^(1-prim.alpha)))./(K_t.^(1-prim.alpha)).-prim.delta  #Initialize underlying functions that will allow for adjustment of w,r
    w_t = ((1-prim.alpha).*(K_t.^(prim.alpha)))./(L_t.^(prim.alpha))
    b_t = (prim.theta.*w_t.*L_t)./sum(F_t[:,:,prim.age_retire:prim.N,1])
    theta = [repeat([0.11], t_iter); repeat([0.0], T-t_iter)]
    res2 = new_primitives(K_t, K_potential, L_t,sav_func_t,val_func_t,lab_func_t, pol_func_t,F_t,pop_normal,Assets,w_t,r_t,b_t,theta) ##Places various objects under res2.
    res2


end
# res2 = Initialize2()




function bellman_t(prim::Primitives,res::Results,res2::new_primitives,T)
     #Function to iterate backwards over ages to determine labor supply, asset holdings, as a function of producitivity, and age and future period asset holding
    #Now we have augmented the function with the transition periods T
    @unpack age_retire, gamma, sigma, z_prod, birth_distribution = prim
    @unpack markov, alpha, delta, beta, N, Assets, na, nz, age_ef, produc = prim
    @unpack w_t, r_t, b_t,theta = res2
    # res.v_final= zeros(na,nz,N)
    # w = 1.05 #######
    # r = 0.05 #######
    # b = 0.2 #######
    # K = 3.3 #######
    # L = 0.3 #######
    # v_0=1e-10
    val_func_t = zeros(na, nz, N, T)
    pol_func_t = zeros(na, nz, N, T)
    lab_func_t= zeros(na, nz, N, T)
    val_func_t[:, :, :, T] = val_func_0 ##Setting initial function values equal to initial steady state
    pol_func_t[:, :, :, T] = pol_func_0
    lab_func_t[:, :, :, T] = lab_func_0
    for t = (T-1):-1:1  ##For all transition time periods
        for a_index = 1:na, z_index= 1:nz ##Over all asset grid values and all states nz
            c = (1+r_t[t])*Assets[a_index] + b_t[t]
            val_func_t[a_index,z_index, N,t] = utility(prim,res,c,0,N) ##Retired aged group
        end

        for age_index = N-1:-1:age_retire ##Retired individuals
            println("Age: ", age_index)
            #elseif age_index >= 46 && age_index < 66  ##The problem for retired agents
            for a_index = 1:na, z_index=1:nz
                val_func_cand = -Inf #Initial value function
                l = 0
            #for a_index = 1:na, z_index = 1:nz
                budget = (1+r_t[t])*Assets[a_index] + b_t[t]
                for ap_index = 1:na, zp_index = 1:nz  ##We now care about saving. Include state just for consistency
                    c = budget - Assets[ap_index]

                    val =utility(prim, res,c,l,age_index)+ beta*val_func_t[ap_index,zp_index,age_index+1,t+1]
                    if val < val_func_cand
                        break
                    end
                    val_func_cand= val
                    pol_func_t[a_index,z_index,age_index,t] = Assets[ap_index] #Update our saving
                    lab_func_t[a_index,z_index,age_index,t] = l #We dont work

                end
                val_func_t[a_index,z_index, age_index,t] = val_func_cand
            end

        end

                    #end
                    # val_func_t[:,:,age_index] .= val
                    # res.v_final[:,:,age_index] .= val
            #else age_index < 46
        for age_index = age_retire-1:-1:1 ##the problem for workers
            for a_index = 1:na, z_index=1:nz
                val_func_cand = -Inf
            #for a_index = 1:na, z_index = 1:nz
                for ap_index = 1:na
                    l = labor_t(prim,res,res2,age_index,a_index,ap_index,z_index,theta,t)

                    budget = (1+r_t[t])*Assets[a_index] + w_t[t]*(1-theta[t])*produc[age_index,z_index]*l  ##Our state is now very important
                    c = budget - Assets[ap_index]
                    val = utility(prim, res,c,l,age_index) .+ beta.*(val_func_t[ap_index,1,age_index+1,t+1].*markov[z_index,1]+val_func_t[ap_index,2,age_index+1,t+1].*markov[z_index,2]) ##value function now heavily impacted by expectation over future productivity states
                    if val < val_func_cand
                        break
                    end
                    val_func_cand= val
                    pol_func_t[a_index,z_index,age_index,t] = Assets[ap_index] #Update our saving
                    lab_func_t[a_index,z_index,age_index,t] = l ##Now we have to set our
                    ##Add the value function over all asset levels for those in the final period of life.

                end
                # val_func_t[:,:,age_index] .= v_0
                val_func_t[a_index,z_index, age_index,t] = val_func_cand

            end
            # res.v_final[:,:,age_index] .= v_0
        end
    end
    val_func_t, pol_func_t, lab_func_t
end


function Fdist_t(prim::Primitives,res2::new_primitives,res::Results,F_0,F_n,T) ##function to get our distribution transition function


    @unpack na,nz,markov,N,Assets,n = prim
    @unpack F = res
    #@unpack T = res2
    F_next= zeros(na,nz,N,T)
    F_next[:,:,:,1] = F_0
    for t_index = 2:T
        F_next[1, :,1, t_index] = F_next[1, :, 1, 1]  ##Must initialize a youngest age goup for all transition periods or else our population mass will leak.
    end
    for t in 1:T-1 ##fror all time periods
        for age_index = 1:N-1
            for a_index in 1:na
                ap_index = 0
                for z_index in 1:nz
                    if F_next[a_index,z_index,age_index,t] > 0
                    #if F_next[a_index,z_index,age_index,T] > 0
                        ap = res2.pol_func_t[a_index,z_index,age_index,t] ##Working backwards to obtain pol func values for particular a,z, and age values in particular transition values t
                        ap_index = 1* argmin(abs.(ap.-Assets)) ##figure out where pol_func value lies on the asset grid
                        for zp_index in 1:nz
                            F_next[ap_index,zp_index,age_index+1,t+1] += markov[z_index,zp_index]*F_next[a_index,z_index,age_index,t]./(1+n) ##Assign F distribution to that asset value
                        #end
                        end
                    end
                end
            end
        end
    end
    F_next

    # for a_index= 1:na,z_index = 1:nz  ###Multiply each measure of the distribution by the normalized population weights
    #     F_next[a_index,z_index,:,:] = F_next[a_index,z_index,:,:].*pop_normal
    # end
    # F_next


end

# res2.F_t = Fdist_t(prim,res2,res,F_0,F_n)
function labor_t(prim::Primitives,res::Results,res2::new_primitives,age_index,a_index,ap_index,z_index,theta,t)  ##Labor Function over transition periods
    @unpack gamma, produc,Assets = prim
    @unpack w_t,r_t,b_t = res2
    l  = (gamma*(1-theta[t])*produc[age_index,z_index]*w_t[t]-(1-gamma)*((1+r_t[t])*Assets[a_index]-Assets[ap_index]))/((1-theta[t])*w_t[t]*produc[age_index, z_index]) ##Optimal Endogenous labor supply, including depreciation
    if l < 0
        l = 0
    elseif l > 1
        l = 1
    end
    l
end

function agg_K_t(prim::Primitives,res2::new_primitives,res::Results,K_0,T)   ##Aggregate capital values over k periods
    @unpack sav_func_t,F_t =res2
    @unpack N, na, nz, Assets = prim

   K_potential = zeros(T)
   K_potential[1] = K_0
   A_array = zeros(na,nz,N)

        # println("We're on Transition period: ",t)
        #for age_index in 1:N
        # println("We're on Age: ",age_index)
    for age_index in 1:N
        # println("We're on Asset: ",a_index)
        for z_index in 1:nz
            A_array[:, z_index,age_index] = Assets
            #K_potential[t+1] += Assets[a_index]*F_t[a_index,z_index,age_index,t+1]
        end
        # sav_func_t[a_index,z_index,age_index,t]
    end
    println(length(A_array))
    for t in 1:T
        K_potential[t] = sum(F_t[:,:,:,t] .* A_array)
    end

    #K_potential = sum(F_t[:,:,:,t].*A_grid)

    return  K_potential
end

        #end
#         K_potential[t]= sum(F_t[:,:,:,t].*A_grid[:,:,:,t])
#     end
#     K_potential
#
# end
# A = new_equilibrium()


# res2.K_potential = new_equilibrium(prim,res2,res,K_0)
# res2.K_potential = agg_L_t(prim,res2,res,K_0)

function agg_L_t(prim::Primitives,res2::new_primitives,res::Results,L_0,T)  ##Aggrefate labor transition path
    @unpack sav_func_t,F_t,lab_func_t =res2
    @unpack N, na, nz, Assets, produc = prim

   L_potential = zeros(T)
   L_potential[1] = L_0
   e_array = zeros(na,nz,N)

        # println("We're on Transition period: ",t)
        #for age_index in 1:N
        # println("We're on Age: ",age_index)
    for a_index = 1:na
        # println("We're on Asset: ",a_index)
        for z_index in 1:nz
            e_array[a_index, z_index,1:45] = produc[:,z_index]
            #K_potential[t+1] += Assets[a_index]*F_t[a_index,z_index,age_index,t+1]
        end
        # sav_func_t[a_index,z_index,age_index,t]
    end

    println(length(e_array))
    for t in 1:T
        L_potential[t] = sum(F_t[:,:,:,t] .* lab_func_t[:,:,:,t] .* e_array)
    end

    #K_potential = sum(F_t[:,:,:,t].*A_grid)

    return  L_potential
end

function KL_update2(prim::Primitives,res2::new_primitives,K_t1,L_t1, lam,F_0,T) #Marking part of the loop a function
    @unpack alpha, theta, age_retire, N, delta,na,nz = prim
    @unpack K_t, L_t, theta = res2
    # retired_mass = sum(pop_normal[age_retire:N])
    K_tnew = lam .* K_t1 .+ (1 - lam) .* K_t
    L_tnew=  lam .* L_t1 .+ (1 - lam) .* L_t
    F_tnew = zeros(na,nz,N,T)
    F_tnew[:,:,:,1] = F_0
    w_tnew = ((1-alpha).*(K_tnew.^(alpha)))./(L_tnew.^(alpha))
    r_tnew = (alpha.*(L_tnew.^(1 .- alpha)))./(K_tnew.^(1 .- alpha)) .- delta
    b_tnew = (theta.*w_tnew.*L_tnew)./sum(F_tnew[:,:,age_retire:N,1])
    res2.K_t = K_tnew
    res2.L_t = L_tnew
    res2.F_t = F_tnew
    res2.r_t = r_tnew
    res2.w_t = w_tnew
    # retired_mass = sum(pop_normal[age_retire:N])
    #for j in age_retired:N #or can do a sum across those numbers in mu, I'm not sure how this will be stored from stat distribution function
     #   retired_mass += (1/1.011)^j
    #end
    #retired_mass = sum(mu(age_retired:N,:,:)) #sum alternative once I know how mass is stored
    res2.b_t = b_tnew
    return
end


###Now We Compute the Transition Path
function overall_solve(prim::Primitives,res::Results,t,T) #Function to solve all code when policy change occurs in t=1
        @unpack alpha, delta, age_retire, N, na,nz = prim
        T_delta = 20
        # T = 30
        tol = 0.01
        i = 0
        lambda = 0.5
        res2 = Initialize2(T,t)
        while true
            while true
                i += 1
                println("***********************************")
                println("Trial #", i)
                res2.val_func_t, res2.pol_func_t, res2.lab_func_t = bellman_t(prim, res, res2, T) #Solve for value, policy and labor transition functions
                res2.F_t = Fdist_t(prim, res2, res, F_0,F_n,T)##get F transition function
                K_t1 = agg_K_t(prim,res2,res,K_0,T) ##Aggregate labor and capital
                L_t1 = agg_L_t(prim,res2,res,L_0,T)
                display(plot([res2.K_t K_t1 repeat([K_0], T) repeat([K_n], T)], ##Displaying capital path over the course of trials
                         label = ["K Guess" "K Path" "Stationary K w/ SS" "Stationary K w/o SS"],
                         title = "Capital Transition Path", legend = :bottomright))
                diff = maximum(abs.(res2.K_t .- K_t1)./K_t1) + maximum(abs.(res2.L_t .- L_t1)./L_t1) #Check gap between guess and second steady state value
                # @printf("Difference: %0.3f.", float(diff))
                println("")
                if diff > tol
                    KL_update2(prim,res2,K_t1,L_t1,lambda,F_0,T)
                    println("Paths Adjusted")
                else
                    println("Paths Converged") ##Once paths have converged print the necessary graphs
                    display(plot([res2.K_t K_t1 repeat([K_0], T) repeat([K_n], T)],
                             label = ["K Guess" "K Path" "Stationary K w/ SS" "Stationary K w/o SS"],
                             title = "Capital Transition Path", legend = :bottomright))
                    savefig("exercise1_aggregate_capital_path.png")
                    println("Aggregate capital path is saved.")
                    println(" ")
                    plot(collect(1:T), [res2.w_t repeat([w_0], T) repeat([w_n], T)],
                      label = ["Wage path" "Stationary wage w/ SS" "Stationary wage w/o SS"],
                      title = "Wage Transition Path", legend = :bottomright)
                      savefig("exercise1_wage_path.png")
                      println("Wage path is saved.")
                      println(" ")
                     plot(collect(1:1:T), [res2.L_t repeat([L_0], T) repeat([L_n], T)],
                         label = ["Aggregate labor path" "Stationary labor w/ SS" "Stationary labor w/o SS"],
                        title = "Labor Transition Path", legend = :bottomright)
                        savefig("exercise1_aggregate_labor_path.png")
                        println("Aggregate labor path is saved.")
                        println(" ")
                    plot(collect(1:T), [res2.r_t repeat([r_0], T) repeat([r_n], T)],
                         label = ["Interest rate path" "Stationary rate w/ SS" "Stationary rate w/o SS"],
                         title = "Interest Rate Transition Path", legend = :bottomright)
                        savefig("exercise1_interest_rate_path.png")
                        println("Interest rate path is saved.")
                        println(" ")



                    break
                end
            end
            error = abs(res2.K_t[T] - K_n)/K_n
            if error > tol
                T += T_delta
                println("Length Increased")
                res2 = Initialize2(T,t)
            else
                println("Iteration Done")
                break
            end
        end
        res2, T
end

elapse = @elapsed res2, T = overall_solve(prim, res,1, 30) ##Run function


function EV(prim::Primitives, res::Results,res2::new_primitives,lam::Float64)  #3Calculating Equivalent Variation values
    @unpack N, na, nz, alpha,sigma = prim
    ev = zeros(na, nz, N)
    evj = zeros(N)
    ev = (res2.val_func_t[:, :, :, 1] ./ val_func_0).^(1/(lam * (1 - sigma)))
    ev_cf = (val_func_n./ val_func_0).^(1/(lam * (1 - sigma)))
    evj_cf = zeros(N)
    for j=1:N
        evj[j] = sum(ev[:, :, j] .* F_n[:, :, j]) ./ sum(F_n[:, :, j])
        evj_cf[j] = sum(ev_cf[:, :, j] .* F_n[:, :, j]) ./ sum(F_n[:, :, j])
    end
    evj, evj_cf, ev, ev_cf
end
evj, evj_cf, ev, ev_cf = EV(prim,res,res2,0.42)
sum((ev.>=1).*F_0)

plot(collect(20:85), [evj[1:66], evj_cf[1:66], repeat([1], 66)], labels = ["With Transition" "Without Transition" "EV = 1"],
          title = "EV by Age", legend = :bottomleft)
savefig("exercise1_EV.png")
println("Average welfare effect within each age cohorts are saved.")
println(" ")
vote_share = sum((ev.>=1).*F_0)
vote_share_counterfact = sum((ev_cf.>=1).*F_0)
diff_vote_share = vote_share_counterfact - vote_share


function overall_solve2(prim::Primitives,res::Results,t,T)  ##function to solve all code when policy occurs in t=21
        @unpack alpha, delta, age_retire, N, na,nz = prim
        T_delta = 20
        # T = 30
        tol = 0.01
        i = 0
        lambda = 0.5
        res2 = Initialize2(T,t)
        while true
            while true
                i += 1
                println("***********************************")
                println("Trial #", i)
                res2.val_func_t, res2.pol_func_t, res2.lab_func_t = bellman_t(prim, res, res2, T)
                res2.F_t = Fdist_t(prim, res2, res, F_0,F_n,T)
                K_t1 = agg_K_t(prim,res2,res,K_0,T)
                L_t1 = agg_L_t(prim,res2,res,L_0,T)
                display(plot([res2.K_t K_t1 repeat([K_0], T) repeat([K_n], T)],
                         label = ["K Guess" "K Path" "Stationary K w/ SS" "Stationary K w/o SS"],
                         title = "Capital Transition Path", legend = :bottomright))
                diff = maximum(abs.(res2.K_t .- K_t1)./K_t1) + maximum(abs.(res2.L_t .- L_t1)./L_t1)
                # @printf("Difference: %0.3f.", float(diff))
                println("")
                if diff > tol
                    KL_update2(prim,res2,K_t1,L_t1,lambda,F_0,T)
                    println("Paths Adjusted")
                else
                    println("Paths Converged")
                    display(plot([res2.K_t K_t1 repeat([K_0], T) repeat([K_n], T)],
                             label = ["K Guess" "K Path" "Stationary K w/ SS" "Stationary K w/o SS"],
                             title = "Capital Transition Path", legend = :bottomright))
                    savefig("exercise2_aggregate_capital_path.png")
                    println("Aggregate capital path is saved.")
                    println(" ")
                    plot(collect(1:T), [res2.w_t repeat([w_0], T) repeat([w_n], T)],
                      label = ["Wage path" "Stationary wage w/ SS" "Stationary wage w/o SS"],
                      title = "Wage Transition Path", legend = :bottomright)
                      savefig("exercise2_wage_path.png")
                      println("Wage path is saved.")
                      println(" ")
                     plot(collect(1:1:T), [res2.L_t repeat([L_0], T) repeat([L_n], T)],
                         label = ["Aggregate labor path" "Stationary labor w/ SS" "Stationary labor w/o SS"],
                        title = "Labor Transition Path", legend = :bottomright)
                        savefig("exercise2_aggregate_labor_path.png")
                        println("Aggregate labor path is saved.")
                        println(" ")
                    plot(collect(1:T), [res2.r_t repeat([r_0], T) repeat([r_n], T)],
                         label = ["Interest rate path" "Stationary rate w/ SS" "Stationary rate w/o SS"],
                         title = "Interest Rate Transition Path", legend = :bottomright)
                        savefig("exercise2_interest_rate_path.png")
                        println("Interest rate path is saved.")
                        println(" ")



                    break
                end
            end
            error = abs(res2.K_t[T] - K_n)/K_n
            if error > tol
                T += T_delta
                println("Length Increased")
                res2 = Initialize2(T,t)
            else
                println("Iteration Done")
                break
            end
        end
        res2, T
end

elapse = @elapsed res2_2, T_2 = overall_solve2(prim, res,20, 50)  ##

evj2, evj_cf2, ev2, ev_cf2 = EV(prim,res,res2_2,0.42) ##Calculating EV values for policy occuring at t= 21
sum((ev2.>=1).*F_0)

plot(collect(20:85), [evj2[1:66], evj_cf2[1:66], repeat([1], 66)], labels = ["With Transition" "Without Transition" "EV = 1"],
          title = "EV by Age", legend = :bottomleft)
savefig("exercise2_EV.png")
println("Average welfare effect within each age cohorts are saved.")
println(" ")
vote_share2 = sum((ev2.>=1).*F_0)
vote_share_counterfact2 = sum((ev_cf2.>=1).*F_0)
diff_vote_share2 = vote_share_counterfact2 - vote_share2

evj2[60]

