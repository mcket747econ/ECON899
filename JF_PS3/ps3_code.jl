# using Parameters, DataFrames, StatFiles, Random, Distributions, Optim, LinearAlgebra

#cd("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS3")
#cd("/Users/jacobbills/Desktop/Economics/Econ 899/JF 3/")

ind_characteristics = DataFrame(load("./Simulated_type_distribution.dta"))
car_ch = DataFrame(load("./Car_demand_characteristics_spec1.dta"))
iv_specification = DataFrame(load("./Car_demand_iv_spec1.dta"))

varlist = ["price","dpm","hp2wt","size","turbo","trans","Year_1986","Year_1987","Year_1988","Year_1989","Year_1990","Year_1991","Year_1992","Year_1993","Year_1994",
"Year_1995","Year_1996","Year_1997","Year_1998","Year_1999","Year_2000","Year_2001","Year_2002","Year_2003","Year_2004","Year_2005",
"Year_2006","Year_2007","Year_2008","Year_2009","Year_2010","Year_2011","Year_2012","Year_2013","Year_2014","Year_2015","model_class_2",
"model_class_3","model_class_4","model_class_5","cyl_2","cyl_4","cyl_6","cyl_8","drive_2","drive_3","Intercept"]

exo_var = varlist[2:length(varlist)]

iv_list = ["i_import","diffiv_local_0","diffiv_local_1","diffiv_local_2","diffiv_local_3","diffiv_ed_0"]




@with_kw struct Param
    T::Int64 = length(unique(car_ch[:,"Year"]))  ##The number of markets(Years)
    # delta_iia::Array{Float64,1} = Vector(car_ch[:,"delta_iia"]) ##the iia delta values
    simdist::Array{Float64,1} = Vector(ind_characteristics[:,"Var1"]) ##Vector containing the simulated individual characteristics
    mz::Array{Float64,1} = Vector(car_ch[:,"price"]) ##Vector containing prices for all models
    eps1::Float64 = 1 
    eps0::Float64 = 10e-12
    years::Array{Float64,1}= sort(unique(car_ch[:,"Year"]))
    share_ch::Array{Float64,1} = Vector(car_ch[:,"share"])
    x::Array{Float64,2} = Matrix(car_ch[:, varlist])
    IV::Array{Float64,2} = hcat(Matrix(car_ch[:,exo_var]),Matrix(iv_specification[:,iv_list]))








end

@with_kw mutable struct Results  ##Include elements that will change over the course of the iterations
    produc_t = Array{Array{Float64,1},1}()
    delta_iia::Array{Float64,1}
    lambda::Float64 = 0.6 ##move to results since this is mutable 
    #res_norm::Array{Float64,1} 






end 


function Initialize() #Initialize parameter struct and results struct
    P = Param()
    produc_t = Array{Array{Float64,1},1}()
    delta_iia = Vector(car_ch[:,"delta_iia"]) ##the iia delta values
    lambda::Float64 = 0.6
    R = Results(produc_t,delta_iia, lambda)

    return P ,R
end 

P,R = Initialize()



function precompute(carch_ch_mt,years,produc_t)  ##Precompute the productids for each particular year
    produc_1 = 0
    for i = 1:31
        produc_1 = findall(x->x==years[i],carch_ch_mt[:,2])
        push!(produc_t,produc_1)

    end
    return produc_t
end


function individual_char(P::Param, produc_t)  
    Zf = Array{Array{Float64,2},1}() ##Structure that will contain all the prices associated with products in a particular market(Year)
    Z = zeros(P.T,1)
    for i=1:P.T
        for j=1:1
            Z = P.mz[Int64.(produc_t[i])]*P.simdist' 
            push!(Zf,Z)
            #Zf[i,j] = Z

        end
    end
    return Zf
end



function value(lambda,t,Z)   ##Evaluating the value of individual characteristics,and prices u
    mu = Array{Float64,2}
    l = 1 ##Number of nonlinear attributes
    mu = lambda.*Z[t]
    mu_exp = exp.(mu)
    return mu_exp
end

function demand(mu,t,delta,produc_t) ##Creating the Demand Function
    ev = exp.(delta[Int64.(produc_t[t])]).*mu
    #ev = exp.(delta[Int64.(produc_t[t])] .- logsumexp.(delta[Int64.(produc_t[t])])).*mu
    market_share = ev./(1 .+ sum(ev,dims=1))
    #println("market share is ", market_share)
    shat = mean(market_share,dims=2) 
    return market_share, shat

end


function jacobian(ms)

    square = ms*(1.0 .- ms)'
    mat2 = ms*ms'#./100
    sz = size(mat2)
    one = ones(sz[1],sz[2])
    id = Diagonal(ones(sz[1],sz[1]))
    #jacob = mean.(eachcol(square)) .- (1/100).*mat2 + Diagonal(zeros(size(mat2)))  #(1/100).*(1 - I).*mat2
    jacob = (1/100) .* id .* square .- (1/100) .* (one-I) .* mat2 #this should be more right
    return jacob
end





function inverse(P::Param,R::Results,eps0,eps1,share_ch,T_thisrun,shat,delta_iia,mu,produc_t,Z,lam)
    delta_res = zeros(size(P.mz,1),1)
    res_norm::Array{Float64,1} = [0]
    if T_thisrun==1
        t= 1 
        mu = value(lam,t,Z) ##Calci
        rowid = produc_t[t]
        f = 100
        i = 1
        jacob = Matrix{Float64}
        while norm(f) > eps1   ##While loop for when norm is greater than 1
            #jacob = 0
            #println("norm f:", norm(f))
            ms, shat= demand(mu,t,R.delta_iia,produc_t)
            #println(i)
            #println("shat is",shat)
            #jacob = jacobian(ms) # not needed in this section
            # println(P.share_ch[Int64.(rowid)])
            # println(log.(shat))
            # println(log.(P.share_ch[Int64.(rowid)]))
            f = log.(P.share_ch[Int64.(rowid)]) .- log.(shat)
            # print(f)
            R.delta_iia[Int64.(rowid)] = R.delta_iia[Int64.(rowid)] .+ f
            println("This is R Delta", R.delta_iia[Int64.(rowid)])
            #println(norm(f))
            append!(res_norm,norm(f))
            #print("res norm really is ", res_norm)
            i += 1
        end 

        while norm(f) > eps0  ##While loop for when norm is less than 1 but above 10e-12
            #jacobian = 1
            println("norm f:", norm(f))
        #println("this is delt ", delta_iia[Int64.(produc_t[T_thisrun])])  
            ev = (R.delta_iia[Int64.(produc_t[t])]).*mu  
        # println("this is delt ", delta_iia[Int64.(produc_t[T_thisrun])])  
            # println("this is ev ", ev)   
            # println("share is", (ev)./(1 .+ ev))
            # println("market share is", mean((ev)./(1 .+ ev),dims=2))
            ms,shat = demand(mu,t,R.delta_iia,produc_t)
            #  ms = (ev)./(1 .+ ev)
            #  shat = mean((ev)./(1 .+ ev),dims=2)
            # println("market share2 is", ms)
        #  println("this is ms ", ms)   
            # println("shat2 is ", shat)
            f = log.(P.share_ch[Int64.(rowid)]).-log.(shat)
            jacob = jacobian(ms)
            # println(" Why is this Na ", R.delta_iia[Int64.(rowid)])
            # println(" Why is this f Na ", f)
            #println(jacob./shat)
            R.delta_iia[Int64.(rowid)] = R.delta_iia[Int64.(rowid)]  + inv(jacob./shat)*f
            println("This is R Delta", R.delta_iia[Int64.(rowid)])
            #println(norm(f))
            append!(res_norm,norm(f))  
            #println("res norm really really is ", res_norm)    
        end
    return res_norm
    
    else
       for t = 1:T_thisrun
            mu = value(lam,t,Z) ##Calci
            rowid = produc_t[t]
            f = 100
            #res_norm::Array{Float64,1} = [0]
            i = 1
            jacob = Matrix{Float64}
            while norm(f) > eps1   ##While loop for when norm is greater than 1
                #jacob = 0
                println("norm f:", norm(f))
                ms, shat= demand(mu,t,R.delta_iia,produc_t)
                #println(i)
                #println("shat is",shat)
                #jacob = jacobian(ms)
                # println(P.share_ch[Int64.(rowid)])
                # println(log.(shat))
                # println(log.(P.share_ch[Int64.(rowid)]))
                f = log.(P.share_ch[Int64.(rowid)]) .- log.(shat)
                # print(f)
                R.delta_iia[Int64.(rowid)] = R.delta_iia[Int64.(rowid)] .+ f
                #println("This is R Delta", R.delta_iia[Int64.(rowid)])
                #println(norm(f))
                append!(res_norm,norm(f))
                #print("res norm really is ", res_norm)
                i += 1
            end 

            while norm(f) > eps0  ##While loop for when norm is less than 1 but above 10e-12
                #jacobian = 1
                println("norm f:", norm(f))
            #println("this is delt ", delta_iia[Int64.(produc_t[T_thisrun])])  
                ev = (R.delta_iia[Int64.(produc_t[t])]).*mu  
            # println("this is delt ", delta_iia[Int64.(produc_t[T_thisrun])])  
                # println("this is ev ", ev)   
                # println("share is", (ev)./(1 .+ ev))
                # println("market share is", mean((ev)./(1 .+ ev),dims=2))
                ms,shat = demand(mu,t,R.delta_iia,produc_t)
                #  ms = (ev)./(1 .+ ev)
                #  shat = mean((ev)./(1 .+ ev),dims=2)
                # println("market share2 is", ms)
            #  println("this is ms ", ms)   
                # println("shat2 is ", shat)
                f = log.(P.share_ch[Int64.(rowid)]).-log.(shat)
                jacob = jacobian(ms)
                # println(" Why is this Na ", R.delta_iia[Int64.(rowid)])
                # println(" Why is this f Na ", f)
                #println(jacob./shat)
                R.delta_iia[Int64.(rowid)] = R.delta_iia[Int64.(rowid)]  + inv(jacob./shat)*f
                #println("This is R Delta", R.delta_iia[Int64.(rowid)])
                #println(norm(f))
                append!(res_norm,norm(f))  
                #println("res norm really really is ", res_norm)    
            end
            delta_res[Int64.(rowid)] = R.delta_iia[Int64.(rowid)]
        end
    end
    #println("res norm is ", res_norm)
    return delta_res
end 






function routine_overall(P::Param,R::Results,t,lam)  ##Overall code function
    P, R = Initialize()
    res_norm::Array{Float64,1} = [0]
    produc_t = Array{Array{Float64,1},1}()
    jacob = Matrix{Float64}
    Zf = Array{Array{Float64,2},1}()
    carch_ch_mt = Matrix(car_ch)
    produc_t = precompute(carch_ch_mt,P.years,produc_t)
    Zf = individual_char(P,produc_t)
    mu =value(lam,t,Zf)
    market_share,shat  = demand(mu,t,R.delta_iia,produc_t)
    #println("shat is ", shat)
    jacob = jacobian(market_share)
    #print(jacob)
    ans = inverse(P,R,P.eps0,P.eps1, P.share_ch,t,shat,R.delta_iia,mu,produc_t,Zf,lam)
   #print(res_norm_init)
   #println("res norm init ", res_norm_init)
    return ans  ##return vector of norms
end

function betaIV(x,inst,W,delta) #delta is the matrix of shats 
    xT = transpose(x)
    iT = transpose(inst)
    xi = xT*inst 
    ix = iT*x
    iD = iT*delta 
    b = inv(xi*W*ix)*xi*W*iD
    return b
end

function resid(P::Param,beta,delta) #delta is the matrix of shats 
    rho = zeros(length(P.share_ch),1)
    for i = 1:length(P.share_ch) #this should cover all j,t combinations since delta is precalculated
        e = P.x[i,:]'*beta
        rho[i] = delta[i] - e[1]
    end
    return rho 
end


function gmmobj(P::Param,R::Results,W,l,t::Int64=31)
    delta = routine_overall(P,R,t,l)
    iT = transpose(P.IV)
    #W = inv(iT*inst), it's better to define W outside the equation, so I can use this for both steps 
    inner = P.IV*W*iT
    b = betaIV(P.x,P.IV,W,delta)
    rho = resid(P,b,delta)
    f = transpose(rho)*inner*rho
    f = f[1]
    return f 
end

function gridsearch(P::Param,R::Results,start,en,n,t::Int64=31)
    L = range(start,en, length = n)
    res = zeros(n,1)
    for l =1:n
        W = inv(transpose(P.IV)*P.IV)
        res[l] = gmmobj(P,R,W,L[l],t) #use the chosen lambda, get the gmm for it
    end
    idx = argmin(res)
    choice = L[idx]
    return res, choice, idx, L #returns the full grid   
end


##last thing to do is the two step

function two_step(P::Param,R::Results,lam_init,t::Int64=31)
    W1 = inv(transpose(P.IV)*P.IV)
    f1(l) = gmmobj(P,R,W1,l,31) #first step 
    opt1 = optimize(f1, [lam_init], BFGS())
    l2 = Optim.minimizer(opt1)
    #get delta with l2
    delta2 = routine_overall(P,R,t,l2)
    b = betaIV(P.x,P.IV,W1,delta2) #might be the wrong weight?
    r2 = resid(P,b,delta2)
    nI = P.IV .* r2 #I assume the intent of the dot in the formula was to signal elementwise   
    W2 = inv(nI'*nI)
    f2(l) = gmmobj(P,R,W2,l,31) #define new objective function with new weights  
    opt2 = optimize(f2, l2, BFGS()) #minimize using those weights 
    l3 = Optim.minimizer(opt2)
    return l2, l3 
end


