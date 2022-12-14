#General Notes for this problem set
##There is no PS1.do on canvas. This is okay becuase estimating a logit is easy, but still a bit annoyin
##Main problem doesn't seem too hard

#Todo

##Solve ML problem

###Use Newton's method

##Solve ML problem again with built in optimizers. 

using StatFiles, DataFrames, LinearAlgebra 

#cd("/Users/jacobbills/Desktop/Economics/Econ 899/JF 1/")

#data = DataFrame(load("Mortgage_performance_data.dta"))

#y = data[!,20]
#y = Array(y)
#ind = ["i_large_loan", "i_medium_loan" , "rate_spread" , "i_refinance", "age_r", "cltv", "dti", "cu" , "first_mort_r", "score_0", "score_1", "i_FHA", "i_open_year2", "i_open_year3", "i_open_year4", "i_open_year5"]
#x = data[!,ind]
#insertcols!(x, 17, :const => 1) #add a constant
#x = Matrix(x)
#β::Array{Float64,1}
#β = zeros(size(x,2)) #I counted sixteen variables in the problem set, plus a constant
#β[17] = -1 #if I understand the problem right

##Question 1 

function LL(beta,y,X) #takes in an input matrix of betas, seems to work 
    N::Int64 = size(y,1) #number of observations
    li::Array{Float64,1} = zeros(N) #initialize output vector
    for i in 1:N
        v = convert(Vector{Float64},X[i,:]) #typing reduces memory allocation and makes this run a lot faster (important for question 4)
        li[i] = -log(1+exp(v'*beta))+y[i]*v'*beta #formula for log-Likelihood
    end
    l::Float64 = sum(li) #log-Likelihood is a single number 
    return l
end

function sx(beta,X) #just a formula I will use a lot
    s::Float64 = 1/(1+exp(-X'*beta))
    return s
end

function score_LL(beta,y,X) #score is first derivatives with respect to beta
    N,K = size(X) #number of observations  
    sl = zeros(K,N) #initialize output array 
    for i in 1:N
        v= convert(Vector{Float64},X[i,:])
        s::Float64 = y[i]-sx(beta,v)  
        sl[:,i] = s*v' # general formula for score 
    end
    score::Array{Float64,1} = zeros(K)
    for k in 1:K
        score[k] = sum(sl[k,:]) #ultimately, want a vector of sums
    end
    return score 
end


function hessian(beta,X) #hessian is second derivatives wrt beta
    N,K = size(X) #getting sizes 
    Hl = zeros(K,K,N) #begin with 3 dimensions
    H = zeros(K,K) #final result is a KxK matrix formed from the sums of each element 
    for i in 1:N
        v = convert(Vector{Float64},X[i,:]) 
        s = sx(beta,v)
        Hl[:,:,i] = v*v'.*s.*(1-s) #most of the formula 
        H += Hl[:,:,i] #getting rid of the number of observations dimension 
    end 
    H = -H #hessians are negative
    return H
end

function evaluate_b(beta,y,X) #runs the above in one function 
    ll = LL(beta,y,X)
    score = score_LL(beta,y,X)
    H = hessian(beta,X)
    return ll, score, H 
end

function newton(beta, y, X, tol::Float64 = 1e-9)
    b0 = beta 
    flag = 0 #convergence flag
    bk::Array{Float64,1} = zeros(size(X,2))
    while flag == 0 #iterate until converged 
        ll, score, H = evaluate_b(b0,y,X) #we only need the score and hessian 
        bk = b0 - H^(-1)*score #formula 
        if norm(bk-b0) < tol 
            b0 = bk 
            flag = 1 #update convergence flag 
        else 
            b0 = bk #No need for a step, though that might make this even faster 
        end 
    end
    ll = LL(b0,y,X) ##care about knowing the maximized log-likelihood 
    return b0, ll #ultimately interested in getting the coefficients and the maximized ll 
end

function q2(beta,y,X)
    N = size(y,1)
    g = zeros(17)
    h = zeros(17,17)
    
    for i in 1:N #get the numerical first and second derivatives evaluated for all observations
        v = convert(Vector{Float64},X[i,:])
        f(b::Vector) = -log(1+exp(v'*b))+y[i]*v'*b #create a function to evaluate at β 
        g[:] += ForwardDiff.gradient(f,beta) #calculates all the partial first derivatives
        h[:,:] += ForwardDiff.hessian(f,beta) #calculates all the partial second derivatives 
    end    
    return g, h 
end