using StatFiles, DelimitedFiles
using Random, Distributions, DataFrames, Parameters, Optim, Plots, LinearAlgebra

#### GHK
# cd("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/JF_PS2") 
##Bring in data
data = DataFrame(load("Mortgage_performance_data1.dta"))


y = data[!,"duration"]
##Set y as our outcome vector


y = Array(y)

ind = ["score_0", "rate_spread", "i_large_loan", "i_medium_loan", "i_refinance", "age_r", "cltv", "dti", "cu", "first_mort_r", "i_FHA", "i_open_year2", "i_open_year3", "i_open_year4", "i_open_year5"]
x = data[!,ind] ##
z = data[!,["score_0", "score_1", "score_2"]] #Set our z vector
insertcols!(x, 16, :const => 1) #add a constant
x = Matrix(x)
z = Matrix(z)
beta = zeros(size(x,2)) #I counted sixteen variables in the problem set, plus a constant
N= size(x,1)



@with_kw struct parameters

    N::Int64 = 16355
    a0::Int64 = 0
    a1::Int64 = -1
    a2::Int64 =-1
    b::Int64 = 0
    g::Float64 =0.3
    rho::Float64 = 0.5
    









end


function Initialize()
    P = parameters()
    return P
end
Random.seed!(1234)

 #setting our seed

#initializing array that will hold our likelihood function


function ll_GHK(P,x)   ##Main GHK Function
    @unpack a0, a1,a2,b,g,rho = P


    sim_matrix = zeros(100,2)  ##Initialize structures we will fill later
    like = zeros(P.N,3,4)
    uni_dist = Uniform(0,1) #Call relevant distributions
    distrib = Normal(0,1)


    

    for i = 1:P.N   ##Looping through all mortgages
        for i =1:P.N
            for t=1:3
               # println(-a0-x[i,:]'*b-z[i,:]'*g)
                like[i,t,1] = cdf(distrib,((-a0-sum(x[i,:]'*b)-sum(z[i,:]'*g))/(1-rho))) #Filling our likelihood function for T=1 
            end
        
        
        end




        for s = 1:100 ##over 100 simulations
            sim_matrix[s,1] = rand(uni_dist)*(1-like[i,1,1])+like[i,1,1] ##Draw epsilon_i0. We truncate by using the complement of the likelihood function for T=1


        end
        


        ###For t= 2

        for t = 1:3
            for s = 1:100
                like[i,t,2] += (1/100) * cdf(distrib,(((-a1-sum(x[i,:]'*b)-sum(z[i,:]'*g)-rho*sim_matrix[s,1]))))
                #Fill our likelihood function for T=2 


            end

        end


        ###For T = 3
        #First we calculate some etas because we will use the autocorrelation structure for epsilons moving forward 
        for t = 1:3
            for s = 1:100
                eta = rand(uni_dist)*(1-like[i,1,2]+like[i,1,2])
                sim_matrix[s,2]= rho*sim_matrix[s,1]+ eta  ##Our new set of epsilons
            end
        end


        for t=1:3
            for s=1:100
                like[i,t,3] += (1/100)*cdf(distrib, ((-a2-sum(x[i,:]'*b)-sum(z[i,:]'*g)-rho*sim_matrix[s,2])))#We fill our likeliohood function for T=3 
            end

        end


        ###For T =4 



        for t=1:3
            for s=1:100
                like[i,t,4] += (1/100)*(1-cdf(distrib, ((-a2-sum(x[i,:]'*b)-sum(z[i,:]'*g)-rho*sim_matrix[s,2])))) #Fill likelihood for T=4
            end

        end
        
        if i % 1000 == 0
            println("We're ", (i/16355)*100, "% Finished")
        end 


    end
 
    return like 


end




function log_like(P,x)  ##Now we solve for the log likelihood 
    like = ll_GHK(P,x) ##Run our  ghk function
    like_step = zeros(P.N,3)
    for i = 1:P.N 
        like_step[i,:] = like[i,:,Int64(y[i])] ##We get the likelihoods in each of the 3 time periods
    end


    total = 0

    like_total = zeros(P.N)
    for i = 1:P.N ##We add all of our likelihoods across each time index together
        if Int64(y[i]) == 1
            like_total[i] = like_step[i,1] 
        elseif Int64(y[i]) == 2
            like_total[i] = like_step[i,2]
        else Int64(y[i]) > 2
            like_total[i] = like_step[i,3]
        end
    end

    for i = 1:P.N 
        total += log(like_total[i]) ##Taking the log, we get the log likelihood for our model. 
    end

    return like, total 

end

P = Initialize()  ##Initialize parameter struct

like,total = log_like(P,x) #Run our overall function and obtain our likelihood array and log likelihood

###Choice Probabilities
#T=1 
mean(like[:,1,1])
#T=2 
mean(like[:,1,2])
#T=3
mean(like[:,1,3])
#T=4
mean(like[:,1,4])


##### Accept-Reject Portion
# cd("/Users/jacobbills/Desktop/Economics/Econ 899/JF 2/")

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

##### Quadrature Portion
# cd("C:\\Users\\Rafeh\\Documents\\GitHub\\ECON899\\PS2_JF")
df = DataFrame(load("Mortgage_performance_data.dta"))
# names(df) |> println

##declare matrices
Y::Array{Float64,2} =  df[:,["i_close_0","i_close_1","i_close_2"]]  |> Matrix

X::Array{Float64,2} = df[:,["score_0","rate_spread","i_large_loan", "i_medium_loan","i_refinance","age_r", "cltv", "dti","cu", "first_mort_r",
 "i_FHA","i_open_year2","i_open_year3","i_open_year4","i_open_year5"]] |>Matrix
X = hcat(ones(size(X)[1]),X)

Z::Array{Float64,2} = df[:,["score_0", "score_1", "score_2"]] |>Matrix
 
## Parameters 
β = zeros(size(X)[2]) 
α = 0,-1,-1
γ = [0.3,0.3,0.3]
ρ = 0.5
σ = 1/abs(1-ρ)


## My understanding
T = zeros(size(Y,1))
T[Y[:,1].==1,1] .=1
T[(Y[:,1].==0) .& (Y[:,2].==1),1] .=2
T[((Y[:,1].==0) .& (Y[:,2].==0)) .& (Y[:,3].==1),1] .=3
T[((Y[:,1].==0) .& (Y[:,2].==0)) .& (Y[:,3].==0),1] .=4


##Likelihood funcitons for t∈{2,3,4}
##1-D - T=2
function lfunct2(x,val,ρ=ρ,σ=σ)
  cdf.(Normal(),(-val-ρ*x)).*
  pdf.(Normal(),(x/σ)/σ)
end 

##2-D - T=3
function lfunct3(x,y,val,ρ=ρ,σ=σ)
  cdf.(Normal(),(-val-ρ*x)).*
  pdf.(Normal(),(x- ρ*y)).*
  pdf.(Normal(),(y/σ)/σ)
end 

##1-D - T=4
function lfunct4(x,y,val,ρ=ρ,σ=σ)
  cdf.(Normal(),(val-ρ*x)).*
  pdf.(Normal(),(x- ρ*y)).*
  pdf.(Normal(),(y/σ)/σ)
end 

##transformation functions
function rho(u,b)
  log(u)+b
end
function rhoprime(u)
  1/u
end

##1-D - T=2
function sumlhood2(x,xmult,w,α,val0)
  sum = 0
  for j in 1:size(x,1)
    sum+=w[j].*lfunct2(x[j]+α[1]+val0,α[2]+val0).*xmult[j]
  end
  sum
end

## T=3 (2-D)
function sumlhood3(x,xmult,w,α,val0)
  sum = 0
  for j in 1:size(x,1)
    # println(rho(x[j,1],val1),", ",rho(x[j,2],val2),", ",rhoprime(x[j,1]),", ",rhoprime(x[j,2])
    sum+=w[j].*lfunct3(x[j,1]+α[2]+val0,x[j,2]+α[1]+val0,α[3]+val0).*xmult[j]
  end
  sum
end

##T=4 (2-D)
function sumlhood4(x,xmult,w,α,val0)
  sum = 0
  for j in 1:size(x,1)
    sum+=w[j].*lfunct4(x[j,1]+α[2]+val0,x[j,2]+α[1]+val0,α[3]+val0).*xmult[j]
  end
  sum
end

k=1
x1 = readdlm("X"*string(k)*".csv",',', Float64)
w1= readdlm("w"*string(k)*".csv",',', Float64)
k=2
x2 = readdlm("X"*string(k)*".csv",',', Float64)
w2 = readdlm("w"*string(k)*".csv",',', Float64)


## Preallocate whatever I can for speed
rhomat_x1 = log.(x1)
x1mult = rhoprime.(x1)

rhoprimemat = rhoprime.(x2)
rhomat_x2 = log.(x2)
x2mult = zeros(size(x2,1))
for k in 1:size(x2mult,1)
  x2mult[k] = rhoprimemat[k,1].*rhoprimemat[k,2]
end


alllhood = function(β,α=α)
  PrT = zeros(size(Y,1))
  val0 = X*β+Z*γ
  for i in 1:size(X,1)
    if T[i]==1
      PrT[i] = cdf.(Normal(),(-α[1]-val0[i])./σ)
    elseif T[i]==2
      PrT[i] = sumlhood2(rhomat_x1,x1mult,w1,α,val0[i])
    elseif T[i]==3
      PrT[i] = sumlhood3(rhomat_x2,x2mult,w2,α,val0[i])
    elseif T[i]==4
      PrT[i] = sumlhood4(rhomat_x2,x2mult,w2,α,val0[i])
    end
  end 
  PrT
end


@time PrT = alllhood(β)
sumlhood = function(β)
  sum(alllhood(β))
end

@time res = optimize(sumlhood, zeros(length(β)),method=LBFGS())

d = DataFrame(T=T,PrT=PrT)
d = hcat(d,DataFrame(Y,:auto))
##
gd = groupby(d, :T)
c = combine(gd, :PrT => mean)
c[!,:norm] =c[!,:PrT_mean]./sum(c[!,:PrT_mean])
c
gd = groupby(DataFrame(T=T,indic=ones(length(T))), :T)
c = combine(gd, :indic => sum)
c[!,:indic] = c[!,:indic_sum]./sum(c[!,:indic_sum])
c |>print

##minimum values
writedlm( "vals_beta.csv",Optim.minimizer(res) , ',')