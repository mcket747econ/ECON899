#General Notes for this problem set
##There is no PS1.do on canvas. This is okay becuase estimating a logit is easy, but still a bit annoyin
##Main problem doesn't seem too hard

#Todo

##Solve ML problem

###Use Newton's method

##Solve ML problem again with built in optimizers. 

using StatFiles, DataFrames, LinearAlgebra 

cd("/Users/jacobbills/Desktop/Economics/Econ 899/JF 1/")

data = DataFrame(load("Mortgage_performance_data.dta"))

y = data[!,20]
y = Array(y)
ind = ["i_large_loan", "i_medium_loan" , "rate_spread" , "i_refinance", "age_r", "cltv", "dti", "cu" , "first_mort_r", "score_0", "score_1", "i_FHA", "i_open_year2", "i_open_year3", "i_open_year4", "i_open_year5"]
x = data[!,ind]
insertcols!(x, 17, :const => 1) 
x = Matrix(x)
β = zeros(size(x,2)) #I counted sixteen variables in the problem set, plus a constant
β[17] = -1 #if I understand the problem right

##Question 1 
function LL(beta,y,X) #takes in an input matrix of betas, seems to work 
    N = size(y,1) #number of observations
    li = zeros(N) #initialize output vector
    for i in 1:N 
        li[i] = -log(1+exp(X[i,:]'*beta))+y[i]*X[i,:]'*beta #formula for log-Likelihood
    end
    l = sum(li) #log-Likelihood is a single number 
    return l
end 

function sx(beta,X) #just a formula I will use a lot 
    1/(1+exp(-X'*beta))
end

function score_LL(beta,y,X) #score is first derivatives with respect to beta
    N,K = size(X) #number of observations  
    sl = zeros(K,N) #initialize output array 
    for i in 1:N
        s = y[i]-sx(beta,X[i,:])  
        sl[:,i] = s*X[i,:]' # general formula 
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
        s = sx(beta,X[i,:])
        Hl[:,:,i] = X[i,:]*X[i,:]'.*s.*(1-s) #most of the formula 
        H += Hl[:,:,i] #getting rid of the number of observations dimension 
    end 
    H = -H
    return H
end

function evaluate_b(beta,y,X) #runs the above in one function 
    ll = LL(beta,y,X)
    score = score_LL(beta,y,X)
    H = hessian(beta,X)
    return ll, score, H 
end

function newton(beta, y, X, tol::Float64 = 1e-9, step::Float64 = .5)
    b0 = beta 
    flag = 0 #convergence flag
    while flag == 0 #iterate until converged 
        ll, score, H = evaluate_b(b0,y,X)
        bk = b0 - H^(-1)*score
        if norm(bk-b0) < tol
            b0 = bk 
            flag = 1 #update convergence flag 
        else 
            b0 = bk #add in a step maybe? 
        end 
    end
    ll = LL(b0,y,X)
    return b0, ll #might want the log-likelihood along with the rest of this
end
