Original_agrid=1; 

@with_kw struct Params
    
    beta::Float64=0.8;
    A::Float64=0.6;

    
    alpha::Float64=0.64; 
    cf::Float64=12; 
    ce::Float64=40; 

    
    tau::Float64=0; 

    p::Float64=1; 
    w::Float64=1; 

    
    rho::Float64=0.93; 
    sigma_logz::Float64=sqrt(0.53); 
    sigma_epsilon::Float64=sqrt((1-rho)*((sigma_logz)^2));
    a::Float64=0.078; 
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
