nes (159 sloc)  5.03 KB

#Original_agrid=1;
using Optim, Plots, Parameters, Distributions, Random, DataFrames

@with_kw struct Params

    beta::Float64=0.8;
    A::Float64=0.6;
    γE::Float64=0.5772156649;

    theta::Float64=0.64;
    cf::Float64=10;
    ce::Float64=5;
    s_grid::Array{Float64,1} = [3.98e-4, 3.58, 6.82, 12.18, 18.79]
    n_s::Int64 = length(s_grid)
    emp_levels::Array{Float64,1} = [1.3e-9, 10, 60, 300, 1000]
    F_transition::Array{Float64,2} = [0.6598 0.2600 0.0416 0.0331 0.0055;
                                      0.1997 0.7201 0.0420 0.0326 0.0056;
                                      0.2000 0.2000 0.5555 0.0344 0.0101;
                                      0.2000 0.2000 0.2502 0.3397 0.0101;
                                      0.2000 0.2000 0.2500 0.3400 0.0100]

    v_s_entrant::Array{Float64,1} =[0.37,0.4631,0.1102,0.0504,0.0063]
    tau::Float64=0;

    p_star::Float64=0.738;
    w::Float64=1;
    lambda::Float64 =.99;


    rho::Float64=0.93;
    sigma_logz::Float64=sqrt(0.53);
    sigma_epsilon::Float64=sqrt((1-rho)*((sigma_logz)^2));
    a::Float64=0.078;
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
end

function Initialize()
    P = Params()
    val_func_in = zeros(P.n_s)
    val_func_out = zeros(P.n_s)
    pf_n_func = zeros(P.n_s)
    pf_prof = zeros(P.n_s)
    pf_entry_x = zeros(P.n_s)
    p_star = 0.738
    mu_0 = ones(P.n_s)/P.n_s
    R = Results(val_func_in,val_func_out,pf_n_func,pf_prof,p_star,pf_entry_x,mu_0)
    return P, R
end
P,R = Initialize()
function n_star(s,theta,p_star,s_grid)
    n = ((1/theta)*p_star*s_grid[s])^((1/theta)-1)
    return n
end
function profit(P::Params,R::Results,i)
    @unpack p_star,pf_n_func, val_func = R
    @unpack s_grid, cf, theta, n_s = P
    #prof_0 = 0
    #for i = 1:n_s
        #n = ((1/theta)*p_star*s_grid[i])^((1/theta)-1)
    prof_1 = p_star*s_grid[i]*n_star(i,theta,p_star,s_grid)^(theta) - n_star(i,theta,p_star,s_grid) - p_star*cf
        #opt = optimize(prof_1,0,100)
        #pf_n[i] = n
        #prof = -opt.minimum()
        #pf_prof[i] = prof
        # if prof> prof_0
        #     prof_0 = prof
        #     pf_prof[i] = prof
        # end

    return prof_1
end
prof = profit(P,R,1)

function Utility(p::Float64,P::Params,R::Results,α::Float64=1)
    @unpack γE= P
    c= 10
    γE/α + (1/α)*(c+log(exp(α*VFI(p,R)-c)))
end

function VFI(P::Params,R::Results)
    @unpack n_s, F_transition, beta= P
    @unpack val_func, pf_entry_x = R
    val_func = zeros(n_s)
    pf_entry_x = zeros(n_s)
    val_func_in = zeros(n_s)
    val_func_out = zeros(n_s)
    for i = 1:n_s
        for i_p = 1:n_s
            val_func_in[i] = profit(P,R,i) + beta*sum(val_func[:].*F_transition[i,:])
            val_func_out[i] += profit(P,R,i)
            if val_func_in[i] > val_func_out[i]
                val_func[i] = val_func_in[i]
                pf_entry_x[i] = 0
            else
                val_func[i] = val_func_out[i]
                pf_entry_x[i] = 1
            end
        end
    end
    return val_func, pf_entry_x
        #Then we want to solve the static labor problem
end

R.val_func, pf_entry_x = VFI(P,R)



function Entval(R::Results,P::Params)
    @unpack n_s = P
    @unpack val_func =  R
    for i = 1:n_s
        W += val_func[i]*v_s_entrant[i]
    end
    return W
end


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
    while convergence_mass = 0
            StatDist(P,R,m0)
        Ld,Ls = LBMC(mu,m,p1)
        if abs(Ld-Ls)>tol
            if Ld>Ls
                m_new = m0 - .01*m0
            elseif Ld <Ls
                m_new = m0 + .01*m0
            end
        elseif abs(Ld-Ls)<tol
            convergemce_mass = 1
        end
    end
    return p1 ## Add return
end
function StatDist(P::Params,R::Results,m,tol::Float64 = 1e-3)
    @unpack n_s,F_transition,v_s_entrant = P
    @unpack mu_0  = R
    error=100

    mu_p = zeros(n_s)
    while error>tol
        for i in 1:n_s
            for ip in 1:n_s
                mu_p[i] += (1-pf_entry_x[i])*F_transition[i,ip]*mu_0[i]+ m*(1-pf_entry_x[i])*F_transition[i,ip]*v_s_entrant[i]
            end
            mu_p
        end
        error = maximum(abs.(mu_p - mu_0))
        mu_0 = mu_p
    end
end
StatDist(P,R,1)


function LBMC(P::Params, R::Results,m)
    @unpack theta, A = P
    @unpack p_star, s_grid, m_0 = R##Labor market clearning
    L_d = 0
    for i = 1:n_s
        L_d += n_star(i,theta,p_star,s_grid)*mu_0[i] +m*(n_star(i,theta,p_star,s_grid)*pf_entry_x[i])
        profit_temp += profit(i)*mu_0[i] + m*(profit(i)*pf_entry_x[i])
    end
    L_s = 1/A - profit_temp
    return L_d, L_s

end

# function Tauchen(mew,sigmasq,rho,znum,q, dshift::Float64=1)
#     sigma=sqrt(sigmasq);
#     zstar=mew/(1-rho);
#     sigmaz=sigma/sqrt(1-rho^2);

#     z=zstar*ones(znum,1) + LinRange(-q*sigmaz,q*sigmaz,znum)';
#     omega=z(2)-z(1);

#     zi=z*ones(1,znum);

#     zj=dshift*ones(znum,znum)+ones(znum,1)*z';

#     P_part1=normcdf(zj+omega/2-rho*zi,mew,sigma);
#     P_part2=normcdf(zj-omega/2-rho*zi,mew,sigma);

#     P=P_part1-P_part2;
#     P(:,1)=P_part1(:,1);
#     P(:,znum)=1-P_part2(:,znum);

#     states=z;
#     transmatrix=P;

#     return states, transmatrix
# end

