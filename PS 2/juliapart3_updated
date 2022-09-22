##This file is my imcomplete notes on how to do part 3

##welfare functions

#need to figure out WFB

#This function should work when we get WFB. I think the notes tell us how to get it, it's for the complete markets
c = 1*(.94) + .5*(1-.4)
α = 1.5
WFB= ((c^(1-α))-1)/(1-α)/(1-.9932)
function lambda(prim::Primitives,res::Results)
    @unpack a_grid, β, δ, α, na, ns, s_grid, markov = prim
    frac = 1/((1-α)*(1-β))
    for s in 1:ns
        for a in 1:na
            lam[a,s] = ((WFB + frac)/(val_func + frac))^(1/(1-α))-1
        end
    end
    return lam
end

#following function needs adjustment but is basically just plugging in formulas
function welfare(prim::Primitives, res::Results, lambda)
    WINC = sum(mu0*val_func)
    WG = sum(lambda*mu0)
    return [WINC,WG]
end

#following function should work fine but needs a working lambda function
function welfare_vote(dist::Distribution)
    vote = zeros(na,ns)
    for i in 1:ns
        for j in 1:na
            if lam[i,j] >= 1
                vote[i,j] = dist[i,j]
            else
                vote[i,j] = 0
            end
        end
    end
    final_count = sum(vote)
    return final_count
end

xyz = lambda(prim,res)
