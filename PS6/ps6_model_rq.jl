#Original_agrid=1;
using Optim, Plots, Parameters, Distributions, Random, DataFrames

@with_kw struct Params
    β::Float64=0.8;
    A::Float64=0.6;     ##Heejin used as 1/200 ?
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
    p_star::Float64=0.738;
    w::Float64=1;
    lambda::Float64 =.99;
    rho::Float64=0.93;
    sigma_logz::Float64=sqrt(0.53);
    sigma_epsilon::Float64=sqrt((1-rho)*((sigma_logz)^2));
    a::Float64=0.078;
    gamma_e::Float64=0.5772156649
     #Distribution of Entrants

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
    p_star = 0.738
    mu_0 = ones(P.n_s)/P.n_s
    R = Results(val_func_in,val_func_out,pf_n_func,pf_prof,p_star,pf_entry_x,mu_0,1,N_d)
    return P, R
end
P,R = Initialize()
function n_star(s,θ,p_star)
    n = ((1/θ)*p_star*s)^((1/θ)-1)
    return n
end
function profit(P::Params,R::Results,s)
    @unpack p_star,pf_n_func, val_func = R
    @unpack cf, θ, n_s = P
    nopt=n_star(s,θ,p_star)
    prof_1 = p_star*s*nopt^(θ) - nopt - p_star*cf
    return prof_1
end
prof = profit(P,R,1)
function Utility(P::Params,R::Results,val_func_0,val_func_1,α::Float64=1)  ##utility
    @unpack γE= P
    @unpack val_func= R
    c= (val_func_0+val_func_1)/2
    γE/α + (1/α)*sum(c+log(exp(α*val_func_0-c) + exp(α*val_func_1-c))), exp(α * (val_func_1 -c))/(exp(α * (val_func_0 - c)) + exp.(α * (val_func_1 - c)))
end

function VFI(P::Params,R::Results)
    @unpack n_s, F_transition, β= P
    val_func = zeros(n_s)
    pf_entry_x = zeros(n_s)
    val_func_in = zeros(n_s)
    val_func_out = zeros(n_s)
    for i = 1:n_s
        val_func_in[i] = profit(P,R,i) + β*sum(val_func[:].*F_transition[i,:])
        val_func_out[i] += profit(P,R,i)
        if val_func_in[i] > val_func_out[i]
            val_func[i] = val_func_in[i]
            pf_entry_x[i] = 0
        else
            val_func[i] = val_func_out[i]
            pf_entry_x[i] = 1
        end
    end
    return val_func, pf_entry_x
        #Then we want to solve the static labor problem
end

function VFI_shocks(P::Params,R::Results,α)
    @unpack n_s, F_transition, β= P
    @unpack pf_prof = R
    val_func = zeros(n_s)
    pf_entry_x = zeros(n_s)
    val_func_1 = zeros(n_s)
    val_func_0 = zeros(n_s)
    for i = 1:n_s
        val_func_0[i] = profit(P,R,i) + β*sum(val_func[:].*F_transition[i,:])
        val_func_1[i] += profit(P,R,i)
        val_func[i],pf_entry_x[i] = Utility(P,R,val_func_0[i],val_func_1[i],α)
    end
    return val_func, pf_entry_x
        #Then we want to solve the static labor problem
end

function VFI_Star(P::Params, R::Results)
    @unpack ns, F, v_s_entrant = P
    @unpack pf_entry_x, val_func, mu, M = R

    muprime = zeros(ns)
    for is = 1:ns
        for isp = 1:ns
            muprime[sp_index] += (1 - pf_entry_x[is]) * mu[is] * 
            F[is, isp] + (1 - pf_entry_x[is]) *
            F[is, isp] * M * v_s_entrant[is]
        end
    end
    muprime
end

# function VFI(P::Params,R::Results)
#     @unpack n_s, F_transition, β= P
#     @unpack val_func, pf_entry_x = R
#     for is = 1:n_s
#         wint = β*sum(val_func[:].*F_transition[is,:])
#         if wint > 0
#             val_func[i] = profit(P,R,i)+wint
#             pf_entry_x[i] = 0
#         else
#             val_func[i] = profit(P,R,i)
#             pf_entry_x[i] = 1
#         end
#     end
#     return val_func, pf_entry_x
# end
R.val_func, pf_entry_x = VFI(P,R)

function Entval(R::Results,P::Params)
    @unpack n_s, v_s_entrant = P
    @unpack val_func =  R
    W=0
    for i = 1:n_s
        W += val_func[i]*v_s_entrant[i]
    end
    return W
end

function StatDist(R::Results,P::Params)
    @unpack pf_entry_x = R
    @unpack F_transition, n_s = P
    for is in 1:n_s
        mu_p = sum((1-pf_entry_x(s))*F_transition[s,s']*mu(s;m))+ m*sum((1-pf_entry_x(s))*F_transition[s,sp]*v_s_entrant)
    end
    mu_p
end


# function LMC(mu,M) ##Labor market clearing
#     @unpack n_s,v_s,s_grid,A = P
#     @unpack mu_0,p_star,M = R
#     L_d::Float64 = 0.0
#     Π::Float64=0.0
#     for is =  1:n_s
#         n_opt = n_star(s_grid[is],p_star)
#         L_d+=n_opt*mu[is] + M*n_opt*v_s[is]
#         Π = profit(P,R)*mu[is] + m*prof(P,R)*v_s[is]
#     end
#     L_s = 1/A - Π
#     return L_d,L_s
# end

function LMC(P::Params,R::Results) ##Labor market clearing
    @unpack n_s,v_s,s_grid,A = P
    @unpack N_d,mu_0,p_star,M = R
    L_d::Float64 = 0.0
    Π::Float64=0.0
    for is =  1:n_s
        n_opt = N_d[is]
        L_d+=n_opt*mu[is] + M*n_opt*v_s[is]
        Π = profit(P,R,s_grid[is])*mu[is] #+ M*prof(P,R)*v_s[is]
    end
    L_s = 1/A - Π
    return L_d,L_s
end


function distr(P::Params, R::Results, tol::Float64 =-1e3)
    @unpack n_s = P
    @unpack pf_entry_x,M = R
    mu = ones(n_s)/n_s
    mu_next = ones(n_s)/n_s
    Wint = 0
    for sp_index=1:n_s
        for s_index = 1:n_s
            Wint+= (1-pf_entry_x[s_index])*F[s_index,sp_index]*mu(s_index)
        end
    end
    mu_next
end


function dist_iterate(prim::Params, res::Results, tol::Float64 = 1e-3)
    error = 100
    n = 0
    while error>tol
        n+=1
        mu_next = VFI_Star(prim,res)
        error = maximum(abs.(mu_next-res.mu))#I assume this is the maximum of the array
        res.mu = mu_next
        println("Iteration number: ", n)
        println("Absolute difference: ", error)
    end
    print("Mu Distribution Converged in ", n, " Iterations.")

end

function solve_firm_prob(R::Results,P::Params,α,tol = 1e-3)
    @unpack n_s, s_grid = P

    n = 1
    err = 100.0

    N_d = zeros(n_s)
    pf_prof = zeros(n_s)
    for is = 1:n_s
        # print(profit(P, R, s_grid[is]))
        N_d[is] = profit(P, R, s_grid[is])
        N_d[is] = n_star(s_grid[is],P.θ,R.p_star)
        pf_prof[is] = profit(P, R, s_grid[is])
    end
    R.N_d = N_d
    R.pf_prof = pf_prof
    while true
        if α == 0.0
            val_func_next, pf_entry_x = VFI(P, R)
        else
            val_func_next, pf_entry_x = VFI_shocks(P, R, α)
        end
        err = maximum(abs.(R.val_func .- val_func_next))
        # println("\n***** ", n, "th iteration *****")
        # @printf("Absolute difference: %0.4f.\n", float(err))
        # println("***************************")
        if err < tol
            break
        end
        R.pf_entry_x = pf_entry_x
        R.val_func = val_func_next
        n += 1

        end
end

function stationary_price(P::Params, R::Results, α::Float64)
    @unpack val_func = R
    @unpack cf, ce,lambda,p_star = P

    p0 = p_star
    converge = 0
    while converge == 0
        R.val_func,pf_entry_x = VFI(P,R)
        W = Entval(R,P)
        if abs(W - p0*ce) > tol
            if W > p0*ce
                p1 = p0 + ((1-p0)/2)
            else
                p1 = p0 - ((1-p0)/2)
            end
        else
            converge = 1
        end
    end
    m0 = m_init
    convergence_mass = 0
    while convergence_mass == 0
        mu = StatDist(R,P)
        Ld,Ls = LMC(mu,m,p1)
        if abs(Ld-Ls)>tol
            if Ld>Ls
                m0-=m_step
            else
                m0+=m_step
            end
        else
            convergence_mass=1
        end
    end
end
##
function solve_HR(P::Params,R::Results,tol=1e-3)
    @unpack val_func = R
    @unpack cf, ce,lambda,p_star = P

    p0 = p_star
    convergence = 0
    while converge == 0
        R.val_func,pf_entry_x = VFI(P,R)
        W = Entval(R,P)
        if abs(W - p0*ce) > tol
            if W > p0*ce
                p1 = p0 + ((1-p0)/2)
            else
                p1 = p0 - ((1-p0)/2)
            end
        else
            convergence = 1
        end
    end
    m0 = m_init
    convergence_mass = 0
    while convergence_mass == 0
        mu = StatDist(R,P)
        Ld,Ls = LMC(mu,m,p1)
        if abs(Ld-Ls)>tol
            if Ld>Ls
                m0-=m_step
            else
                m0+=m_step
            end
        else
            convergence_mass=1
        end
    end
    return p1,mu,m0 ## Add return
end

function solve_HR2(P::Params,R::Results,tol=1e-3)
    @unpack val_func = R
    @unpack cf, ce,lambda,p_star = P

    p0 = p_star
    convergence = 0
    while converge == 0
        R.val_func,R.pf_entry_x = VFI(P,R)
        W = Entval(R,P)
        if abs(W - p0*ce) > tol
            if W > p0*ce
                p1 = p0 + ((1-p0)/2)
            else
                p1 = p0 - ((1-p0)/2)
            end
        else
            convergence = 1
        end
    end
    R.m = 1
    Ld,Ls = LMC(mu,m,p1)
        if abs(Ld-Ls)>tol
            if Ld>Ls
                m0-=m_step
            else
                m0+=m_step
            end
        else
            convergence_mass=1
        end
    return p1,mu,m0 ## Add return
end


function EVS()
    u_0 = zeros(n_s)
    for i =1:n_s
        val_func_in = profit(P,R,i) + β *sum(u_0[:].*F_transition[i,:])
        val_funcOut = profit(P,R,i)
    end
    U_1[i] = (gamma_e/alpha_1) +(1/alpha_1)*log((e^(alpha_1*val_func_in[i])+e^(alpha_1*val_func_out[i])))
end