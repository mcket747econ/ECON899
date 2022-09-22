# using Pkg; Pkg.add(["Plots","Parameters"])
using Plots, Parameters
@with_kw struct Primitives
    β::Float64 = 0.9932 #Discount Rate
    b::Float64 = 0.5 #Unemployment Benefits
    α::Float64 = 1.5
    price::Array{Float64,1} = collect(range(0,length=1000,stop=1))
    S::Array{Float64,1} = [1, 0.5]
    Π::Array{Float64,2} = [0.97 .03; 0.5 0.5]
    A:: Array{Float64,1}= collect(range(-2,length=1000,stop=5))
    ne::Int64 = length(S)
    na::Int64 = length(A)


end

mutable struct Results
    val_func:: Array{Float64,2}
    pol_func:: Array{Float64,2}
    mu:: Array{Float64,2}
    q:: Float64
end


function Initialize()
    prim = Primitives()
    val_func = zeros(prim.na,prim.ne)
    pol_func = zeros(prim.na,prim.ne)
    mu = ones(prim.na, prim.ne)/(prim.na*prim.ne)## The division is in order to create a uniform distribution
    q = (1 + prim.β)/2
    res = Results(val_func, pol_func, mu, q)
    # V_iterate(prim,res)
    prim, res



end
function V_iterate(prim::Primitives, res::Results, tol::Float64 = 1e-3)
    error = 100
    n = 0
    while error>tol
        n+=1
        v_next = Bellman(prim,res)
        error = maximum(abs.((v_next-res.val_func))) #I assume this is the maximum of the array
        res.val_func = v_next
        println("Iteration number: ", n)
        println("Absolute difference: ", error)
    end
    print("Value Function Converged in ", n, " Iterations.")

end

function Bellman(prim::Primitives,res::Results)
    @unpack  Π, β, na, ne, price, α, b, A, S = prim
    v_next = zeros(na, ne)

    for s_index=1:ne
        ap_index_start = 1
        for a_index = 1:na
            v_0 = -1e10
            budget = S[s_index] + A[a_index]
            for ap_index=ap_index_start:na
                c = budget - A[ap_index]*res.q
                if c > 0
                    v = ((c^(1-α))-1)/(1-α) + β*sum(res.val_func[ap_index,:].*Π[s_index,:])
                    if v > v_0
                        v_0 = v
                        res.pol_func[a_index,s_index] = A[ap_index]
                        ap_index_start = ap_index
                    end
                end
            end
            v_next[a_index,s_index] = v_0
        end
    end
    v_next
end

function distribution(prim::Primitives, res::Results, tol::Float64 =-1e3)
    @unpack Π, β, na, ne, price, α, b, A, S = prim
    mu_next = zeros(na, ne)
    # mu_0 = ((A[0]-(-2))/(7))
    for a_index=1:na,s_index=1:ne
        ap = res.pol_func[a_index,s_index]
        ap_index = argmin(abs.(ap .- A))
        #s_p = res.pol_func(a_index,:)
        for sp_index = 1:ne
            mu_next[ap_index,sp_index] += res.mu[a_index, s_index] * Π[s_index, sp_index]
        end

    end
    mu_next
end


function dist_iterate(prim::Primitives, res::Results, tol::Float64 = 1e-3)
    error = 100
    n = 0
    while error>tol
        n+=1
        mu_next = distribution(prim,res)
        error = maximum(abs.(mu_next-res.mu))#I assume this is the maximum of the array
        res.mu = mu_next
        println("Iteration number: ", n)
        println("Absolute difference: ", error)
    end
    print("Mu Distribution Converged in ", n, " Iterations.")

end

function price_update(prim::Primitives, res::Results, tol::Float64 = 1e-3)
    aggregate_a = sum(res.mu[:,1] .* prim.A ) + sum(res.mu[:,2] .* prim.A ) #if it is positive: oversaving, negative: overborrowing.
    if abs(aggregate_a) > tol
        if aggregate_a > 0
            res.q = res.q + (1 - res.q)/2 * abs(aggregate_a)
        elseif aggregate_a < 0
            res.q = res.q - (res.q-prim.β)/2 * abs(aggregate_a)
        end
    end
end

prim, res = Initialize()


#V_iterate(prim,res)
#dist_iterate(prim,res)
#price_update(prim, res)

function overall_iterate(prim::Primitives, res::Results, tol::Float64 = 1e-3) #doing 1e-4 gets stuck in a non working loop
    balance = 100
    x = 0
    while abs(balance) > tol
        x+=1
        V_iterate(prim,res)
        dist_iterate(prim,res)
        balance = aggregate_a = sum(res.mu[:,1] .* prim.A ) + sum(res.mu[:,2] .* prim.A )
        price_update(prim,res)
        println("Iteration number: ", x)
        println("Absolute difference: ", balance)
    end


end

##gini index stuff

#need a function to create wealth distribution, below should work as mu is an (N,2) array
function wealth(dist,prim::Primitives) #for turning mu into a wealth distribution
    @unpack na, A, S = prim
    #first create wealth grids for each state
    wgrid_e = zeros(na,2)
    wgrid_u = zeros(na,2)
    wgrid_e[:,1] = A .+ S[1]
    wgrid_u[:,1] = A .+ S[2]
    #attach probabilities to each wealth spot
    wgrid_e[:,2] = dist[:,1]
    wgrid_u[:,2] = dist[:,2]
    #combine the two arrays into one array
    wealth_dist = sortslices(vcat(wgrid_e,wgrid_u),dims = 1) #okay to do this only because no overlap in wealth between the groups
    return wealth_dist
end

#goal of the above is to turn the two state mu distribution into a single vector of assets plus earnings

function gini(dist) #no need to sort, the distribution is sorted already. Special thanks to Wikipedia and Stack Exchange
    #first column of distribution is wealth, second is probability
    nrows = size(dist,1)
    S_n = zeros(nrows) #creating a vector for Sn (see wikipedia)
    for i in 1:nrows #iterate over the different wealths
        for j in 1:i #but only up to a certain amount each time
        S_n[i] += dist[j,2]*dist[j,1] #summing the weighted wealth over wealth up to that point
        end
    end

    num = 0
    for i in 2:nrows
        num += dist[i,2]*(S_n[i-1]+S_n[i])
    end

    gini_in = 1 - num/(S_n[length(S_n)])
    return gini_in
end

##Lorenz curve

#function plot_Lorenz(dist) #input an ordered distribution
 #   n_d = size(dist,1)
  #  lorenz = zeros(n_d,2)
   # lorenz[:,1] = collect(range(0, length = n_d, stop = 1)) #goes from 0 to 1
   # lorenz[1,2] = dist[1,1]*dist[1,2] #create the first value
   # total = 0
   # count = 0
   # for i in 2:n_d #loop over rest of curve
    #     total = dist[i,1]*dist[i,2] + lorenz[i-1,2] #amount of wealth at gridpoint plus a cumulative summation
     #    lorenz[i,2] = total
      #   count = i
       # if total > .99 #because the bounds aren't working right with all the zeros
        #    break
       # else
        #    continue
       # end
    #end
    #Plots.plot(lorenz[1:count,1],lorenz[1:count,2], title="Lorenz curve") #basic Lorenz plot
    #Plots.plot!(lorenz[1:count,1],lorenz[1:count,1]) #45 degree line
#end


overall_iterate(prim,res)

Plots.plot(prim.A,res.pol_func, title = "Decision Rules", label = ["Employed" "Unemployed"])
Plots.plot!(prim.A,prim.A) #45 degree line on above plot

Plots.plot(prim.A, res.mu, title="Asset Distribution", label = ["Employed" "Unemployed"])

wf = wealth(res.mu,prim)
print("Gini coefficient: ", gini(wf))
#plot_Lorenz(wf) #Plots a Lorenz curve.


#sum(res.mu[:,])
sum(res.mu[:,1])
sum(res.mu[:,2])


#end
c = 1*(0.9433962264150745) + .5*(1-0.9433962264150745)

WFB= ((c^(1-α))-1)/(1-α)/(1-.9932)
function lambda(prim::Primitives,res::Results,WFB)
    @unpack A, β, α, na, ne, S, Π = prim
    @unpack val_func, mu = res
    frac = 1/((1-α)*(1-β))
    lam = zeros(na,ne)
    for s in 1:ne
        for a in 1:na
            lam[a,s] = ((WFB + frac)/(val_func[a,s] + frac))^(1/(1-α))-1
        end
    end
    return lam
end

#following function needs adjustment but is basically just plugging in formulas
function welfare(prim::Primitives, res::Results, lambda)
    @unpack val_func, mu = res
    WINC = sum(mu.*val_func)
    WG = sum(lambda.*mu)
    return [WINC,WG]
end

#following function should work fine but needs a working lambda function
function welfare_vote(prim::Primitives,res::Results,lam)
    @unpack ne, na = prim
    @unpack mu = res
    vote = zeros(na,ne)
    for i in 1:na
        for j in 1:ne
            if lam[i,j] >= 0
                vote[i,j] = mu[i,j]
            else
                vote[i,j] = 0
            end
        end
    end
    final_count = sum(vote)
    return final_count
end

xyz = lambda(prim,res,WFB)
sum(res)|

Plots.plot(prim.A, xyz, title = "Consumption Equivalence", label = ["Employed" "Unemployed"])
Plots.savefig("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/Problem Set 1/02Growthzipjulia/02_policy_Functions.png")

lambda
welf = welfare(prim,res,xyz)
print(welf)
print(xyz)
vote = welfare_vote(prim,res,xyz)
