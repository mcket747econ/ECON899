using Plots, Parameters
using Pkg;Pkg.add("UnPack")
@with_kw struct Primitives
    beta::Float64 = .9932
    alpha::Float64= 1.5
    S::Array{Float64,1} = [1, 0.5]
    markov::Array{Float64,2} = [.97 0.03;0.5 0.5]
    Assets:: Array{Float64,1} = collect(range(-2,length=1000,stop=5))
    na::Integer = length(Assets)
    ns::Integer = length(S)

end



@with_kw mutable struct Results
    val_func::Array{Float64,2}
    pol_func::Array{Float64,2}
    q_price::Float64
end

function Initialize()
    prim = Primitives()
    val_func = zeros(prim.na,prim.ns)
    pol_func = zeros(prim.na,prim.ns)
    q_price = (1+ prim.beta)/2 #Choosing a starting price that is between 1 and the our discount rate
    res= Results(val_func,pol_func,q_price) #I'm still not entirely sure what is happening here
    prim, res #We return these structures so that we can access them in the rest of our programs

end

function V_iterate(prim::Primitives,res::Results,tol=1e-3)
    error = 100 #We begin with a high error
    n = 0
    while error > tol
        n += 1
        v_next = Bellman(prim,res)
        error = maximum(abs.(v_next - res.val_func))
        res.val_func = v_next

        println("Iteration Number: ", n)
        println("Absolute Difference:", error)
    end
    println("Function Converged In: ", n, "Iterations.")
end


function Bellman(prim::Primitives,res::Results)
    @unpack  markov, beta, na, ns, alpha, Assets, S = prim
    v_next = zeros(prim.na,prim.ns)
    #I'm not sure where to put the initial guess
    for a_index = 1:prim.na,s_index = 1:prim.ns #First we iterate over assets and states
        budget = prim.Assets[a_index] + prim.S[s_index]
        ap_index_start=1
        v_0=-1e10 #The goal here is for an initial asset level, determine the optimal amount of next period assets
        for ap_index = ap_index_start:prim.na,sp_index = 1:prim.ns ##Now3 we begin iterating over possible future states
            c = prim.Assets[a_index] + prim.S[s_index] - prim.Assets[ap_index]*res.q_price  #We now come to the aprim indices because our consumption i impacted by tomorow's choice
            if c > 0
                v = (c^(1-prim.alpha)-1)/(1-prim.alpha) +prim.beta*sum(res.val_func[ap_index,:].*prim.markov[s_index,:])
                if v > v_0
                    v_0 = v
                    res.pol_func[a_index,s_index] = prim.Assets[ap_index]
                    ap_index_start = ap_index
                end
            end
        end
        v_next[a_index,s_index] = v_0 #for some reason we need to stay inside the loop
    end
    v_next
end

prim, res = Initialize()
V_iterate(prim,res)
Plots.plot(prim.A,res.val_func,title = 'Self-Test Graph',label=['Good State' 'Bad State'])
Plots.plot(prim.A,res.pol_func,title= "Self Test Policy Graph", label=['Good State' 'Bad State'])








end
