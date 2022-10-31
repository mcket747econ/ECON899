#General Notes for this problem set
##There is no PS1.do on canvas. This is okay becuase estimating a logit is easy, but still a bit annoyin
##Main problem doesn't seem too hard

#Todo

##Write Three functions (all conditional on \beta)

###Log-Likelihood

###Score of Log-Likelihood

###Hessian (annoying but okay)

##Evaluate at two different betas {-1,0}

##Solve ML problem

###Use Newton's method

##Solve ML problem again with built in optimizers. 

using StatFiles, DataFrames

cd("/Users/jacobbills/Desktop/Economics/Econ 899/JF 1/")

data = DataFrame(load("Mortgage_performance_data.dta"))

y = data[!,20]
y = Array(y)
ind = ["i_large_loan", "i_medium_loan" , "rate_spread" , "i_refinance", "age_r", "cltv", "dti", "cu" , "first_mort_r", "score_0", "score_1", "i_FHA", "i_open_year2", "i_open_year3", "i_open_year4", "i_open_year5"]
x = data[!,ind]
insertcols!(x, 17, :const => 1) 
x = Matrix(x)
β = zeros(size(x)) #I counted sixteen variables in the problem set, plus a constant
β[17] = -1 #if I understand the problem right

##Question 1 
function LL(beta,y,X) #takes in an input matrix of betas, seems to work 
    N = size(y,1) #number of observations
    li = zeros(N)
    for i in 1:N
        li[i] = -log(1+exp(X[i,:]'*beta)+y[i]*X[i,:]'*beta)
    end
    l = sum(li)
    return l
end 

function sx(beta,X)
    1/(1+exp(-X'*beta))
end

function score_LL(beta,y,X) #score is first derivatives with respect to beta
    N,K = size(X)
    sl = zeros(K,N)
    for i in 1:N
        s = y[i]-sx(beta,X[i,:])
        sl[:,i] = s*X[i,:]'
    end
    score::Array{Float64,1} = zeros(K)
    for k in 1:K
        score[k] = sum(sl[k,:])
    end
    return score 
end


function hessian(beta,X) #hessian is second derivatives wrt beta
    N,K = size(X)
    Hl = zeros(K,K,N)
    H = zeros(K,K)
    for i in 1:N
        s = sx(beta,X[i,:])
        Hl[:,:,i] = X[i,:]*X[i,:]'.*s.*(1-s)
        H += Hl[:,:,i]
    end 
    return H
end

function evaluate_b(beta,y,X)
    ll = LL(beta,y,X)
    score = score_LL(beta,y,X)
    H = hessian(beta,X)
    return ll, score, H 
end