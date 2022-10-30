#Original_agrid=1;
using Optim, Plots, Parameters, Distributions, Random, DataFrames

@with_kw struct Params

    beta::Float64=0.8;
    A::Float64=0.6;


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



function Entval(R::Results
    @unpack val_func, R
  for i = 1:n_s
      W += val_func[i]*v_s_entrant[i]
  end
  return W
end


function solve_HR(P::Params,R::Results,tol=)
    @unpack val_func = R
    @unpack cf, ce,lambda = P

    p0 = p_star
    convergence = 0
    while converge == 0
        R.val_func,pf_entry_x = VFI()
        W = Entval(R)
        if abs(W - p0*ce) > tol
            if W > p0*ce
                p1 = p0 + ((1-p0)/2)
            else
                p1 = p0 - ((1-p0)/2)
        else
            convergence = 1

    end

    m0 = m_init
    convergence_mass = -
    while convergence_mass = 0







end

function StatDist()
    mu_p = sum((1-pf_entry_x(s))*F_transition[s,s']*mu(s;m))+ m*sum((1-pf_entry_x(s))*F_transition[s,sp]*v_s_entrant)




end

function Tauchen(mew,sigmasq,rho,znum,q, tauchenoptions, dshift)

    sigma=sqrt(sigmasq);
    zstar=mew/(1-rho);
    sigmaz=sigma/sqrt(1-rho^2);

    z=zstar*ones(znum,1) + linspace(-q*sigmaz,q*sigmaz,znum)';
    omega=z(2)-z(1);

    zi=z*ones(1,znum);

    zj=dshift*ones(znum,znum)+ones(znum,1)*z';

    P_part1=normcdf(zj+omega/2-rho*zi,mew,sigma);
    P_part2=normcdf(zj-omega/2-rho*zi,mew,sigma);

    P=P_part1-P_part2;
    P(:,1)=P_part1(:,1);
    P(:,znum)=1-P_part2(:,znum);

    states=z;
    transmatrix=P;

    return states, transmatrix
end



n_z=20;
Params.q=4;
[z_grid, pi_z]=TauchenMethod(Params.a,Params.sigma_epsilon^2,Params.rho,n_z,Params.q,tauchenoptions);
z_grid=exp(z_grid);

n_a=601;
a_grid=[linspace(0,100,101),((logspace(0,pi,n_a-101)-1)/(pi-1))*(5000-101)+101]';
if ImposeFootnote5==1
    a_grid=[linspace(0,100,101),((logspace(0,pi,n_a-101-1)-1)/(pi-1))*(5000-101)+101,10^6]';
    if Original_agrid==1
        n_a=250;
        a_grid=[exp(linspace(0,log(5001),n_a-1))-1,10^6]';
    end
end
n_d=0;
d_grid=[];




if ImposeFootnote5==1
    vfoptions.ReturnToExitFn=@(a_val, z_val,tau) -tau*a_val*(a_val~=10^6);
end


[V,Policy,ExitPolicy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);




figure(1)
surf(ExitPolicy)




Params.upsilon=[1;zeros(n_a-1,1)]*[ones(1,floor(0.65*n_z)),zeros(1,n_z-floor(0.65*n_z))];
Params.upsilon=Params.upsilon/sum(sum(Params.upsilon));


simoptions
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions,Params,EntryExitParamNames);








if ImposeFootnote5==0
    plot(a_grid, cumsum(sum(StationaryDist.pdf,2)))
else
    temp=sum(StationaryDist.pdf,2);
    temp2=temp(1:end-1); temp2(1)=temp(1)+temp(end);
    plot(a_grid(1:end-1), cumsum(temp2))
end



GEPriceParamNames={'ce','Ne'};

FnsToEvaluateParamNames(1).Names={'alpha'};



FnsToEvaluateFn_1 = @(aprime_val,a_val,z_val,agentmass,alpha) z_val*(aprime_val^alpha);
FnsToEvaluate={FnsToEvaluateFn_1};


AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);



heteroagentoptions.specialgeneqmcondn={0,'entry'};


GeneralEqmEqnParamNames(1).Names={'p','A'};

GeneralEqmEqn_1 = @(AggVars,GEprices,p,A) A/AggVars-p;
GeneralEqmEqnParamNames(2).Names={'p','beta'};
GeneralEqmEqn_Entry = @(EValueFn,GEprices,p,beta) beta*EValueFn-p*GEprices(1);
if ImposeFootnote5==1
    GeneralEqmEqnParamNames(2).Names={'p'};
    GeneralEqmEqn_Entry = @(EValueFn,GEprices,p) EValueFn-p*GEprices(1);
end



GeneralEqmEqns={GeneralEqmEqn_1, GeneralEqmEqn_Entry};


heteroagentoptions.verbose=1;
n_p=0;
disp('Calculating price vector corresponding to the stationary eqm')


[p_eqm,p_eqm_index, GeneralEqmCondition]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions, EntryExitParamNames);

[V,Policy,ExitPolicy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
Params.zeta=1-ExitPolicy;
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions,Params,EntryExitParamNames);

