using StatFiles, DelimitedFiles
using Random, Distributions, DataFrames, Parameters, Optim, Plots, LinearAlgebra
pyplot()

cd("C:\\Users\\Rafeh\\Documents\\GitHub\\ECON899\\PS2_JF")
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

writedlm( "vals_beta.csv",Optim.minimizer(res) , ',')
=