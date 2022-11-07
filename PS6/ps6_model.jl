
using Optim, Plots, Parameters, Distributions, Random, DataFrames

@with_kw struct Params
    β::Float64=0.8;
    A::Float64=1/200;     ##Heejin used as 1/200 ?
    γE::Float64=0.5772156649;

    θ::Float64=0.64;
    cf::Float64=10;
    ce::Float64=5;
    s_grid::Array{Float64,1} = [3.98e-4, 3.58, 6.82, 12.18, 18.79]
    n_s::Int64 = length(s_grid)
    F_transition::Array{Float64,2} = [0.6598 0.2600 0.0416 0.0331 0.0055;
                                      0.1997 0.7201 0.0420 0.0326 0.0056;
                                      0.2000 0.2000 0.5555 0.0344 0.0101;
                                      0.2000 0.2000 0.2502 0.3397 0.0101;
                                      0.2000 0.2000 0.2500 0.3400 0.0100]

    v_s_entrant::Array{Float64,1} =[0.37,0.4631,0.1102,0.0504,0.0063]
    
    emp_levels::Array{Float64,1} = [1.3e-9, 10, 60, 300, 1000]
    tau::Float64=0;
    p_star::Float64=1.0;
    w::Float64=1;
    lambda::Float64 =.99;
    rho::Float64=0.93;
    sigma_logz::Float64=sqrt(0.53);
    sigma_epsilon::Float64=sqrt((1-rho)*((sigma_logz)^2));
    a::Float64=0.078;
    gamma_e::Float64=0.5772156649
end

mutable struct Results
    val_func::Array{Float64,1}
    val_func_test::Array{Float64,1}
    pf_n_func::Array{Float64,1}
    pf_prof::Array{Float64,1}
    p_star::Float64
    pf_entry_x::Array{Float64,1}
    mu_0::Array{Float64,1}
    M::Float64
    N_d::Array{Float64,1}
end

function Initialize()
    P = Params()
    val_func_in = zeros(P.n_s)
    val_func_out = zeros(P.n_s)
    pf_n_func = zeros(P.n_s)
    pf_prof = zeros(P.n_s)
    pf_entry_x = zeros(P.n_s)
    N_d = zeros(P.n_s)
    p_star = 1.0
    mu_0 = ones(P.n_s)     #/P.n_s
    R = Results(val_func_in,val_func_out,pf_n_func,pf_prof,p_star,pf_entry_x,mu_0,5.0,N_d)
    return P, R
end
P,R = Initialize()

function n_star(s,θ,p_star)
    n = (θ*p_star*s)^(1/(1-θ))
    if n < 0
        n=0
    end
    return n
end

function profit(P::Params,R::Results,s)
    @unpack p_star = R
    @unpack cf, θ = P
    nopt=n_star(s,θ,p_star)
    prof_1 = p_star*s*nopt^(θ) - nopt - p_star*cf
    return prof_1
end

function Utility(P::Params,R::Results,val_func_0,val_func_1,α::Float64=1)  ##utility
    @unpack γE= P
    @unpack val_func= R
    c= (val_func_0+val_func_1)/2
    γE/α + (1/α)*(c+log(exp(α*val_func_0-c) + exp(α*val_func_1-c))), 
    exp(α * (val_func_1 -c))/(exp(α * (val_func_0 - c)) + exp(α * (val_func_1 - c)))
end

function VFI(P::Params,R::Results)
    @unpack n_s, F_transition, β= P
    @unpack val_func, pf_prof = R

    pf_entry_x_out = zeros(n_s)
    val_func_in = zeros(n_s)
    val_func_out = zeros(n_s)
    val_func_next = zeros(n_s)
    for i = 1:n_s
        val_func_out[i] = pf_prof[i]
        val_func_in[i] =  pf_prof[i]+ β*sum(val_func[:].*F_transition[i,:])
        if val_func_in[i] >= val_func_out[i]
            val_func_next[i] = val_func_in[i]
            pf_entry_x_out[i] = 0
        else
            val_func_next[i] = val_func_out[i]
            pf_entry_x_out[i] = 1
        end
    end
    return val_func_next, pf_entry_x_out
end

function VFI_shocks(P::Params,R::Results,α)
    @unpack n_s, F_transition, β= P
    @unpack pf_prof, pf_entry_x,val_func = R
    # pf_entry_tmp = zeros(n_s)
    val_func_1 = zeros(n_s)
    val_func_0 = zeros(n_s)
    for i = 1:n_s
        val_func_0[i] =pf_prof[i] + β*sum(val_func[:].*F_transition[i,:])
        val_func_1[i] = pf_prof[i]
        val_func[i],pf_entry_x[i] = Utility(P,R,val_func_0[i],val_func_1[i],α) 
    end
    return val_func, pf_entry_x
        #Then we want to solve the static labor problem
end

function VFI_Star(P::Params, R::Results)
    @unpack n_s, F_transition, v_s_entrant = P
    @unpack pf_entry_x, val_func, mu_0, M = R

    muprime = zeros(n_s)
    for is = 1:n_s
        for isp = 1:n_s
            muprime[isp] += (1 - pf_entry_x[is]) * mu_0[is] * 
            F_transition[is, isp] + (1 - pf_entry_x[is]) *
            F_transition[is, isp] * M * v_s_entrant[is]
        end
    end
    muprime
end


function Entval(R::Results,P::Params)
    @unpack n_s, v_s_entrant = P
    @unpack val_func =  R
    W=0
    for i = 1:n_s
        W += val_func[i]*v_s_entrant[i]
    end
    return W
end

function StatDist(R::Results,P::Params,tol=1e-3)
    @unpack n_s = P
    n=0
    mu_p =zeros(n_s)
    while true
        # mu_p = sum((1-pf_entry_x(s))*F_transition[s,s']*mu(s;m))+ m*sum((1-pf_entry_x(s))*F_transition[s,sp]*v_s_entrant)
        mu_p= VFI_Star(P,R)
        diff = maximum(abs.(mu_p.-R.mu_0))
        if diff<tol 
            break
        end 
        R.mu_0 = mu_p
        n+=1
    end
    return mu_p
end


function LMC(P::Params,R::Results) ##Labor market clearing
    @unpack n_s,v_s_entrant,s_grid,A = P
    @unpack N_d,mu_0,p_star,M,pf_prof = R
    L_d::Float64 = 0.0
    Π::Float64=0.0
    for is =  1:n_s
        n_opt = N_d[is]
        L_d+=n_opt*mu_0[is] + M*n_opt*v_s_entrant[is]
        Π += pf_prof[is]*mu_0[is] #+ M*prof(P,R)*v_s[is]
    end
    L_s = 1/A - Π
    return L_d,L_s
end


function solve_firm_prob_shocks(R::Results,P::Params,α,tol = 1e-3)
    @unpack n_s, s_grid = P

    n = 1
    err = 100.0

    N_d = zeros(n_s)
    pf_prof = zeros(n_s)
    for is = 1:n_s
        # print(profit(P, R, s_grid[is]))
        N_d[is] = n_star(s_grid[is],P.θ,R.p_star)
        pf_prof[is] = profit(P, R, s_grid[is])
    end
    R.N_d = N_d
    R.pf_prof = pf_prof
    while true
        val_func_next, pf_entry_x = VFI_shocks(P, R, α)
        err = maximum(abs.(R.val_func .- val_func_next))
        # println("\n***** ", n, "th iteration *****")
        # @printf("Absolute difference: %0.4f.\n", float(err))
        # println("***************************")
        if err < tol
            R.val_func = val_func_next
            break
        end
        R.pf_entry_x = pf_entry_x
        R.val_func = val_func_next
        n += 1
        end
end



function solve_firm_prob_certain(R::Results,P::Params,tol = 1e-3)
    @unpack n_s, s_grid = P

    n = 1
    err = 100.0

    N_d = zeros(n_s)
    pf_prof = zeros(n_s)
    for is = 1:n_s
        # print(profit(P, R, s_grid[is]))
        N_d[is] = n_star(s_grid[is],P.θ,R.p_star)
        pf_prof[is] = profit(P, R, s_grid[is])
    end
    R.N_d = N_d
    R.pf_prof = pf_prof
    while true
        val_func_next, pf_entry_x = VFI(P, R)
        err = maximum(abs.(R.val_func .- val_func_next))
        # println("\n***** ", n, "th iteration *****")
        # @printf("Absolute difference: %0.4f.\n", float(err))
        # println("***************************")
        if err < tol
            R.val_func = val_func_next
            break
        end
        R.pf_entry_x = pf_entry_x
        R.val_func = val_func_next
        n += 1
    end
end


function stationary_price(P::Params,R::Results,α::Float64,tol=1e-3)
    @unpack v_s_entrant, ce,n_s = P
    @unpack val_func, p_star = R
    n  = 1
    converge=0
    while converge==0
        if α==0.0
            solve_firm_prob_certain(R,P,tol)
        else
            solve_firm_prob_shocks(R,P,α,tol)
        end
        EC = 0.0
        for is in 1:n_s
            EC+= R.val_func[is] *v_s_entrant[is]
        end
        EC= EC/R.p_star - ce
        if abs(EC)<tol 
            break
        end 
        if EC>=0 
            R.p_star=R.p_star* (1-1e-4)
        else
            R.p_star=R.p_star* (1+1e-4)
        end
        if n%200==0
            println("p is ",R.p_star,"    -    n is ",n)
            println("tol is ",tol,"     -    EC is ",EC)
        end
        if n==100000 
            printl("Kill Switch")    ##Kill switch @ 17k
            converge=1
            break
        end
        n+=1
    end
end


function stationary_M(P::Params,R::Results,tol=1e-3)
    @unpack val_func,mu_0 = R
    @unpack cf, ce, lambda, p_star = P
    convergence_mass = 0
    n=0
    # M_next=0.0
    while convergence_mass == 0
        StatDist(R,P)
        Ld,Ls = LMC(P,R)
        if abs(Ld-Ls)>tol
            if Ld>Ls
                R.M-=R.M*1e-4
            else
                R.M+=R.M*1e-4
            end
        else
            convergence_mass=1
            println("Convergence Achieved")
            break
        end
        if n % 200 == 0
            println("************************************")
            println("Absolute value ", float(abs(Ld-Ls)))
            println("=> M is update ", float(R.M))
        end
        if n == 50000 
            println("Reached 1000")
            break
        end
        n+=1
    end
    return R ## Add return
end


## Table to export 
table= function(P::Params,R::Results)
    @unpack n_s,v_s_entrant = P
    mass_stay = sum(R.mu_0 .* (1 .- R.pf_entry_x))
    mass_exit = sum(R.mu_0 .* R.pf_entry_x)
    L_stay = sum(R.N_d .* R.mu_0)
    L_enter = sum(R.N_d * R.M .* v_s_entrant)
    agg_lab = sum(R.N_d .* R.mu_0  + R.N_d * R.M .* v_s_entrant)
    [R.p_star, R.M, mass_stay,mass_exit, agg_lab, L_stay, L_enter, L_enter/agg_lab]
end

