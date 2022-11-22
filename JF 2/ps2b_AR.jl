#General Notes for this problem set
##There is no PS2.do on canvas. This is okay becuase estimating a probit is easy, but still a bit annoyin

using StatFiles, DataFrames, LinearAlgebra,  Distributions, Random 

cd("/Users/jacobbills/Desktop/Economics/Econ 899/JF 2/")

data = DataFrame(load("Mortgage_performance_data.dta"))

y = data[!,"duration"] 
y = Array(y)
ind = ["score_0", "rate_spread", "i_large_loan", "i_medium_loan", "i_refinance", "age_r", "cltv", "dti", "cu", "first_mort_r", "i_FHA", "i_open_year2", "i_open_year3", "i_open_year4", "i_open_year5"]
x = data[!,ind]
z = data[!,["score_0", "score_1", "score_2"]]
insertcols!(x, 16, :const => 1) #add a constant
x = Matrix(x)
z = Matrix(z)
beta::Array{Float64,1} = zeros(size(x,2)) #I counted sixteen variables in the problem set, plus a constant
N::Int64 = size(x,1)


alpha0::Int64 = 0
alpha1::Int64 = -1
alpha2::Int64 = -1
gamma::Array{Float64,1} = [0.3,0.3,0.3]
rho::Float64 = 0.5


##Question 3

function draw_e(N, r, n::Int64=100, seed=12345)
    e::Array{Float64,3} = zeros(N,3,n)
    Random.seed!(seed)
    for x in 1:n 
        e[:,1,x] = rand(Normal(0,1/((1-r)^2)),N)
        e[:,3,x] = rand(Normal(0,1),N)
        e[:,2,x] = r.*e[:,1,x] + e[:,3,x]
    end
    return e 
end

function ar2(e,x,z,b,g,a0,a1,n,r,s)
    R = 0
    count = 0
    for w in 1:n
        if e[w] < a0+ x'*b + z'*g
            C = cdf.(Normal(),(-a1-x'*b-z'*g-r*e[w]))
            p = pdf.(Normal(),(e[w]/s))/s
            R += C*p 
            count += 1
        end
    end
    m = (1/count)*R
    return m 
end

function ar3(e,x,z,b,g,a0,a1,a2,n,r,s)
    R = 0
    count = 0
    for w in 1:n
        if e[2,w] < a1+ x'*b + z'*g && e[1,w] < a0+ x'*b + z'*g
            C = cdf.(Normal(),(-a2-x'*b-z'*g-r*e[2,w]))
            p = pdf.(Normal(),(e[1,w]/s))/s
            f = pdf.(Normal(),(e[2,w]-r*e[1,w]))
            R += C*f*p
            count += 1
        end
    end
    m = (1/count)*R
    return m 
end

function ar4(e,x,z,b,g,a0,a1,a2,n,r,s)
    R = 0
    count = 0
    for w in 1:n
        if e[2,w] < a1+ x'*b + z'*g && e[1,w] < a0+ x'*b + z'*g
            C = cdf.(Normal(),(a2+x'*b+z'*g-r*e[2,w]))
            p = pdf.(Normal(),(e[1,w]/s))/s
            f = pdf.(Normal(),(e[2,w]-r*e[1,w]))
            R += C*f*p
            count += 1
        end
    end
    m = (1/count)*R
    return m 
end


function ll_ar(b,g, x, z, y, N, a0::Int64 = 0, a1::Int64 = -1,a2::Int64 = -1, r::Float64 = 0.5, n::Int64=100)
    e = draw_e(N,r,n) #get simulated error terms
    s = 1/sqrt(1-r)
    P::Array{Float64,1} = zeros(N)
    ll = 0 
    prob = 0
    pseudo::Array{Float64,2} = zeros(N,5)
    for i in 1:N
        pseudo[i,1] = cdf.(Normal(0,1),(-a0-x[i,:]'*b-z[i,:]'*g)/s)
        pseudo[i,2] = ar2(e[i,1,:],x[i,:],z[i,:],b,g,a0,a1,n,r,s)
        pseudo[i,3] = ar3(e[i,:,:],x[i,:],z[i,:],b,g,a0,a1,a2,n,r,s)
        pseudo[i,4] = ar4(e[i,:,:],x[i,:],z[i,:],b,g,a0,a1,a2,n,r,s)
        if y[i] == 1
            P[i] = pseudo[i,1]
            pseudo[i,5] = 1
        elseif y[i] ==2
            P[i] = pseudo[i,2]
            pseudo[i,5] = 2
        elseif y[i] ==3
            P[i] = pseudo[i,3]
            pseudo[i,5] = 3
        else
            P[i] = pseudo[i,4]
            pseudo[i,5] = 4
        end
        prob = -log(max(0,P[i]))
        ll += prob 
    end
    return ll, pseudo 
end

 res, prob = ll_ar(beta,gamma,x,z,y,N)

 println("Average Probability T=1 over all observations: ", mean(prob[:,1]))
 println("Average Probability T=2 over all observations: ", mean(prob[:,2]))
 println("Average Probability T=3 over all observations: ", mean(prob[:,3]))
 println("Average Probability T=4 over all observations: ", mean(prob[:,4]))
