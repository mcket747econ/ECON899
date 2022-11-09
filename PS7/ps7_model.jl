using Random, Distributions, DataFrames, Parameters, Optim, Plots, LinearAlgebra
pyplot()


@with_kw struct parameters
    rho_0::Float64 = 0.5
    sigma_0::Float64 = 1
    x_0::Float64 = 0
    T::Int64 = 200
    l::Int64 = 2
    H_0::Int64 = 1
    H::Int64 = 10
    W = Matrix{Int64}(I,3,3)
    f_rho::Array{Float64,1} = collect(range(0.35,length=100,stop=0.65))
    f_sig::Array{Float64,1} = collect(range(0.8,length=100,stop=1.2))
    grid_l::Int64 = 100







end

@with_kw mutable struct results
    x_t::Array{Float64,1}
    y_t::Array{Float64,2}
    e_0::Array{Float64,2}
    e::Array{Float64,2}
    M_t::Array{Float64,1}
    #M_th::Array{Float64,1}
    j_func::Array{Float64,2}
    rho_hat::Float64
    sig_hat::Float64
     #shocks::Array{Float64,1}






end

#onstruct(results,e_0::Array{Float64,1})
function initialize()
    P = parameters()
    x_t = zeros(P.T)
    y_t = zeros(P.T,P.H)
    e_0 = rand(Normal(0,(P.sigma_0)^2),P.T,P.H_0)
    e = rand(Normal(0,(P.sigma_0)^2),P.T,P.H)
     #shocks = zeros(P.T)
    M_t = zeros(2)
    #M_th = (1/(P.T*P.H))
    j_func = zeros(P.grid_l,P.grid_l)
    rho_hat = 0
    sig_hat = 0
    R = results(x_t,y_t,e_0,e,M_t,j_func,rho_hat,sig_hat)


    return P,R
end

Random.seed!(1234)
P,R = initialize()


function simulate(P,e,f_sig,f_rho)
    #set.seed(123)
    x_t = zeros(P.T,P.H)
    x_t[1,:] = f_rho*P.x_0 .+  f_sig*e[1,:]
    for i in 2:P.T
        x_t[i,:] = f_rho*x_t[i-1,:] .+ f_sig*e[i,:]
    end

    return x_t

end

function simulate_0(P,e)
    #set.seed(123)
    x_t = zeros(P.T,P.H_0)
    x_t[1,:] = P.rho_0*P.x_0 .+ e[1,:]
    for i in 2:P.T
        x_t[i,:] = P.rho_0*x_t[i-1,:] .+ P.sigma_0*e[i,:]
    end

    return x_t

end

###Initial_simulation
function testing()
    xy= simulate_0(P,R.e_0)
    R.x_t = xy[:,1]

    x_tm1 = zeros(200)
    for t=1:P.T
        if t > 1
            x_tm1[t] = R.x_t[t-1]
        else
            x_tm1[t] = 0
        end
    end
    R.M_t = [mean(R.x_t), sum(((R.x_t .- mean(R.x_t)).*((R.x_t .- mean(R.x_t)))))/(P.T*P.H_0), sum(((R.x_t .- mean(R.x_t)).*((x_tm1 .- mean(R.x_t)))))/(P.T*P.H_0)]
end


# function Getmoments(P,f_sig,f_rho,R,T,H)
#     #@unpack y_t = R
#     y_t = zeros(T,H)
#     y_tm1 = zeros(T,H)
#     if H ==1
#         y_t = simulate_0(P,R.e_0)
#     else
#         y_t = simulate(P,R.e,f_sig,f_rho)
#     end
#     m_2 = zeros(2,H)
#     for t=1:T
#         if t > 1
#             y_tm1[t] = y_t[t-1]
#         else
#             y_tm1[t] = 0
#         end
#     end
#
#     for i = 1:H
#         m_2[:,i] = [sum(((y_t .- mean(y_t)).*((y_t .- mean(y_t)))))/(P.T*H), sum(((y_t .- mean(y_t)).*((y_tm1 .- mean(y_t)))))/(P.T*H)]
#     end
#     M_th = mean(m_2,dims=2)
#     return M_th
# end

###GetMoments for 3 moment case
function Getmoments(P,f_sig,f_rho,R,T,H)
    #@unpack y_t = R
    y_t = zeros(T,H)
    y_tm1 = zeros(T,H)
    if H ==1
        y_t = simulate_0(P,R.e_0)
    else
        y_t = simulate(P,R.e,f_sig,f_rho)
    end
    m_2 = zeros(3,H)
    for t=1:T
        if t > 1
            y_tm1[t] = y_t[t-1]
        else
            y_tm1[t] = 0
        end
    end

    for i = 1:H
        m_2[:,i] = [mean(y_t),sum(((y_t .- mean(y_t)).*((y_t .- mean(y_t)))))/(P.T*H), sum(((y_t .- mean(y_t)).*((y_tm1 .- mean(y_t)))))/(P.T*H)]
    end
    M_th = mean(m_2,dims=2)
    return M_th
end





### Get Moments for 2 moment case, mean and variance
# function Getmoments(P,f_sig,f_rho,R,T,H)
#     #@unpack y_t = R
#     y_t = zeros(T,H)
#     if H ==1
#         y_t = simulate_0(P,R.e_0)
#     else
#         y_t = simulate(P,R.e,f_sig,f_rho)
#     end
#     m_2 = zeros(2,H)
    # for t=1:T
    #     if t > 1
    #         y_tm1[t] = y_t[t-1]
    #     else
    #         y_tm1[t] = 0
    #     end
    # end
#     for i = 1:H
#         m_2[:,i] = [mean(y_t), sum(((y_t .- mean(y_t)).*((y_t .- mean(y_t)))))/(P.T*H)]
#     end
#     M_th = mean(m_2,dims=2)
#     return M_th
# end


####Getmoments for two moment case, variance autocovariance

# function Getmoments(P,f_sig,f_rho,R,T,H)
#     #@unpack y_t = R
#     y_t = zeros(T,H)
#     if H ==1
#         y_t = simulate_0(P,R.e_0)
#     else
#         y_t = simulate(P,R.e,f_sig,f_rho)
#     end
#     R.y_t = simulate(P,R.e,f_sig,f_rho)

    # for t=1:T
    #     if t > 1
    #         y_tm1[t] = y_t[t-1]
    #     else
    #         y_tm1[t] = 0
    #     end
    # end
#     m_2 = zeros(2,H)
#     for i = 1:H
#         m_2[:,i] = [sum(((y_t .- mean(y_t)).*((y_t .- mean(y_t)))))/(P.T*P.H),sum(((y_t .- mean(y_t)).*((y_tm1 .- mean(y_t)))))/(P.T*H)]
#     end
#     M_th = mean(m_2,dims=2)
#     return M_th
# end

# function objective
#     M_th = Getmoments()
#     j_func = sum((M_t .- M_th).*W.*(M_t.-M_th))
#
#
#
# end


# function SMM_fm(x,M_t,W,e,P)
#     P.f_rho = x(1)
#     P.f_sigma = x(2)
#     obj = SMM(M_t,W,e,P)
# end
#
# function minobjfunc(W,k)
#
#
#
# end

# function gamma(P,R,j)
#     #@unpack y_t = R
#     y_t = R.y_t
#     M_th = zeros(2,P.H)
#     m_2 = zeros(2,P.H)
#     #y_t = zeros(P.T,P.H)
#     gamm = zeros(2)
#     for i = 1:P.H
#         m_x  = mean(y_t[:,i])
#         m_2[:,i] = [mean(y_t), sum((y_t .- mean(y_t)).*((y_t .- mean(y_t))/P.T))]
#     end
#     M_th = mean(m_2,dims=2)
#     for h = 1:P.H
#         xbar = mean(y_t[:,h])
#         for t = j+1:P.T
#             gamm = gamm + (1/((P.T)*P.H)).*(([y_t[t,h], (y_t[t,h].-xbar)^2] .- M_th).*([y_t[t-j,h],(y_t[t-j,h].-xbar)^2  ] .- M_th))
#         end
#     end
#     return gamm
#
# end



#R.e_0[1][0:5,1]

###underidentified Case
function objective(P,M_t,W,e,R)
    M_th = zeros(P.grid_l)
    j_func = zeros(P.grid_l,P.grid_l)
    for i = 1:P.grid_l
        #P.f_sig = P.f_sig[i]
        for j = 1:P.grid_l
            #f_rho = P.f_rho[j]
            #j_func[i,j],M_th[i,j] = SMM(M_t,W,e,P,P.f_sig[i],P.f_rho[j])
            M_th = Getmoments(P,P.f_sig[i],P.f_rho[j],R,P.T,P.H)
            # print(M_th)
            # print(M_th)
            j_func[i,j] = sum((M_t .- M_th).*W.*(M_t .- M_th))
        end
    end
    return j_func
end
Getmoments(P,P.f_sig[50],P.f_rho[1],R,P.T,10)
R.M_t
    #M_th = zeros(P.grid_l,P.grid_l)
R.j_func = objective(P,R.M_t,P.W,R.e,R)
function min_obj(j_func)
    sig_hat_index = argmin(j_func)[1]
    #println(sig_hat_index)
    rho_hat_index = argmin(j_func)[2]
    #println(rho_hat_index)
    sig_hat = P.f_sig[sig_hat_index]
    rho_hat = P.f_rho[rho_hat_index]

    return sig_hat,rho_hat
end



plot(P.f_rho,P.f_sig,R.j_func,st=:surface,xlabel ="Rho",ylabel="Sigma",zlabel="J Function",camera=(40,2))
savefig("C:/Users/mcket/OneDrive/Documents/Fall 2022/ECON899-Computational/899Code/PS7/j_func5.png")
sig_hat,rho_hat = min_obj(R.j_func)
argmin(R.j_func)
minimum(R.j_func)
function ComputeSE(e,rho_hat,sigma_hat,P,R,H,W)
    eps = 1e-12
    m_2 = Getmoments(P,sigma_hat,rho_hat,R,P.T,H)
    m_rhop = Getmoments(P,sigma_hat,rho_hat .- eps,R,P.T,H)
    # print(m_rhop)
    #println(m_rhop)
    m_sigp = Getmoments(P,sigma_hat.-eps,rho_hat,R,P.T,H)
    delta_rho = (m_rhop .- m_2)./eps
    delta_sig = (m_sigp .- m_2)./eps
    delta_rho1 = delta_rho[[1,2,3]]
    # println("delta_rho1",delta_rho1)
    delta_b = [delta_rho,delta_sig]
    println("deltab ",delta_b)
    delta_b = [delta_b[1],delta_b[2]]
    # println("deltab ",delta_b)

    x = delta_b'*I
    y = x.*delta_b
    # println("y",y)
    y2 = y[[1,2,3]]
    # println(y)
    # println("y2",y2)
    # typeof(y)
    vcov =  inv.(y)
    # print("vcov",vcov)
    #println(vcov)
    return (diag(vcov)).^(0.5)


end

min_obj(R.j_func)
delta_b = [[-0.9655887200921143; -30.63105324940807;;], [-0.13414269695033454; -8.033129716977783;;]]
delta_b[2]
function procedure_smm(P,R,lag)
    m_0 = Getmoments(P,P.sigma_0,P.rho_0,R,P.T,1)
    #println("mo is",m_0)
    R.j_func = objective(P,m_0,P.W,R.e,R)
    rho_hat, sigma_hat = min_obj(R.j_func)
#    println("first rho is ",rho_hat)
#    println("first sigma is ",sigma_hat)
    diag_vcov = ComputeSE(R.e,rho_hat,sigma_hat,P,R,P.H,P.W)
    m_1 = Getmoments(P,sigma_hat,rho_hat,R,P.T,P.H)
    #println(m_1)
    y_t = simulate(P,R.e,sigma_hat,rho_hat)
    S_th = NeweyWest(P,m_1,y_t)
    W_star = inv(S_th)
    #println('W', W_star)
    R.j_func = objective(P,m_1,W_star,R.e,R)
    rho_hat_2, sigma_hat_2 = min_obj(R.j_func)
#    println(rho_hat_2)
#    println("sigma is",sigma_hat_2)
    diag_vcov_2 = ComputeSE(R.e,rho_hat_2,sigma_hat_2,P,R,P.H,W_star)
    return rho_hat_2,sigma_hat_2,diag_vcov_2,diag_vcov

end
argmin(R.j_func)
minimum(R.j_func)*200*(10/11)
procedure_smm(P,R,2)
inv([[939.1937847444417;;] [246.1927507952068;;]; [246.1927507952068;;] [64.54916731293666;;]])
min_obj(R.j_func)
R.j_func
Getmoments(P,P.sigma_0,P.rho_0,R,P.T,1)

[[0.0010647429915351322;;] [0.0040618580229108405;;]; [0.0040618580229108405;;] [0.015492066615700314;;]]
# function NeweyWest(m_0,m_1,lag,P)
#         S_yth = gamma(P,R,0)
#         for j =1:lag
#             gamm = gamma(P,R,0)
#             S_yth = S_yth + (1 - (j/(P.T+1)))*(gamm .+ permutedims(gamm))
#         end
#         S_th = (1+(1/P.H))*S_yth
#
#         return S_th
# end

x = [[6.3211015; -1.5558487],[5.8445315;-1.491324795357]]
y = x'*I
x
w = [0.07981572046047004 0.00047865492685299156; 0.00047865492685299156 0.017884976566461413]
a = x'*I
c = a.*x
xyz = [[939.1937847444417;;] [246.1927507952068;;]; [246.1927507952068;;] [64.54916731293666;;]]
[xyz[[1,2]]
y= x.*w
y*x
function NeweyWest(prim::parameters, m_0,y_t)
    lag_max = 4
    Sy = GammaFunc(prim, m_0,y_t, 0)

    # loop over lags
    for i = 1:lag_max
        gamma_i = GammaFunc(prim, m_0,y_t, i)
        Sy += (gamma_i + gamma_i').*(1-(i/(lag_max + 1)))
    end
    S = (1 + 1/prim.H).*Sy

    return S
end

Getmoments(P,P.sigma_0,P.rho_0,R,P.T,1)
R.j_func
Getmoments(R.e,P,.8,.8,R,P.T,P.H)


###GammaFunc for overidentified case
function GammaFunc(prim::parameters, m_0,y_t, lag::Int64)
    @unpack H, T = prim

    mom_sim = [m_0[1], m_0[2],m_0[3]]
    data_sim = y_t

    # gamma_tot = zeros(length(prim.grid_l),length(prim.grid_l))
    gamma_tot = zeros(3,3)

    for t = (1+lag):T
        for h = 1:H
            # No Lagged
            avg_obs = data_sim[t,h]
            if t > 1
                avg_obs_tm1 = data_sim[t-1,h]
            else
                avg_obs_tm1 = 0
            end
            avg_h = mean(data_sim[:,h])
            var_obs = (avg_obs - avg_h)^2
            auto_cov_obs = (avg_obs - avg_h)*(avg_obs_tm1 - avg_h)

            mom_obs_diff = [avg_obs, var_obs, auto_cov_obs] - mom_sim
            # mom_obs_diff = [avg_obs, var_obs, auto_cov_obs] - mom_sim
            mom_obs_diff = mom_obs_diff

            # Lagged
            avg_lag = data_sim[t-lag,h]
            if t - lag > 1
                avg_lag_tm1 = data_sim[t-lag-1,h]
            else
                avg_lag_tm1 = 0
            end
            avg_h = mean(data_sim[:,h])
            var_lag = (avg_lag - avg_h)^2
            auto_cov_lag = (avg_lag - avg_h)*(avg_lag_tm1 - avg_h)

            mom_lag_diff = [avg_lag, var_lag,auto_cov_lag] - mom_sim
            # mom_lag_diff = [avg_lag, var_lag, auto_cov_lag] - mom_sim
            mom_lag_diff = mom_lag_diff
            #print(mom_lag_diff)
            #print(mom_obs_diff)

            gamma_tot = gamma_tot .+ mom_obs_diff*mom_lag_diff'
        end
    end

    gamma = (1/(T*H)).*gamma_tot

    return gamma
end

####GammaFunc for just identified mean and variance case
function GammaFunc1(prim::parameters, m_0,y_t, lag::Int64)
    @unpack H, T = prim

    mom_sim = [m_0[1], m_0[2]]
    data_sim = y_t

    # gamma_tot = zeros(length(prim.grid_l),length(prim.grid_l))
    gamma_tot = zeros(2,2)

    for t = (1+lag):T
        for h = 1:H
            # No Lagged
            avg_obs = data_sim[t,h]
            if t > 1
                avg_obs_tm1 = data_sim[t-1,h]
            else
                avg_obs_tm1 = 0
            end
            avg_h = mean(data_sim[:,h])
            var_obs = (avg_obs - avg_h)^2
            auto_cov_obs = (avg_obs - avg_h)*(avg_obs_tm1 - avg_h)

            mom_obs_diff = [avg_obs, var_obs] - mom_sim
            # mom_obs_diff = [avg_obs, var_obs, auto_cov_obs] - mom_sim
            mom_obs_diff = mom_obs_diff

            # Lagged
            avg_lag = data_sim[t-lag,h]
            if t - lag > 1
                avg_lag_tm1 = data_sim[t-lag-1,h]
            else
                avg_lag_tm1 = 0
            end
            avg_h = mean(data_sim[:,h])
            var_lag = (avg_lag - avg_h)^2
            auto_cov_lag = (avg_lag - avg_h)*(avg_lag_tm1 - avg_h)

            mom_lag_diff = [avg_lag, var_lag] - mom_sim
            # mom_lag_diff = [avg_lag, var_lag, auto_cov_lag] - mom_sim
            mom_lag_diff = mom_lag_diff
            #print(mom_lag_diff)
            #print(mom_obs_diff)

            gamma_tot = gamma_tot .+ mom_obs_diff*mom_lag_diff'
        end
    end

    gamma = (1/(T*H)).*gamma_tot

    return gamma
end
##gamma func for just idenftied variance autocovariance case
function GammaFunc2(prim::parameters, m_0,y_t, lag::Int64)
    @unpack H, T = prim

    mom_sim = [m_0[2],m_0[3]]
    data_sim = y_t

    # gamma_tot = zeros(length(prim.grid_l),length(prim.grid_l))
    gamma_tot = zeros(3,3)

    for t = (1+lag):T
        for h = 1:H
            # No Lagged
            avg_obs = data_sim[t,h]
            if t > 1
                avg_obs_tm1 = data_sim[t-1,h]
            else
                avg_obs_tm1 = 0
            end
            avg_h = mean(data_sim[:,h])
            var_obs = (avg_obs - avg_h)^2
            auto_cov_obs = (avg_obs - avg_h)*(avg_obs_tm1 - avg_h)

            mom_obs_diff = [ var_obs, auto_cov_obs] - mom_sim
            # mom_obs_diff = [avg_obs, var_obs, auto_cov_obs] - mom_sim
            mom_obs_diff = mom_obs_diff

            # Lagged
            avg_lag = data_sim[t-lag,h]
            if t - lag > 1
                avg_lag_tm1 = data_sim[t-lag-1,h]
            else
                avg_lag_tm1 = 0
            end
            avg_h = mean(data_sim[:,h])
            var_lag = (avg_lag - avg_h)^2
            auto_cov_lag = (avg_lag - avg_h)*(avg_lag_tm1 - avg_h)

            mom_lag_diff = [var_lag,auto_cov_lag] - mom_sim
            # mom_lag_diff = [avg_lag, var_lag, auto_cov_lag] - mom_sim
            mom_lag_diff = mom_lag_diff
            #print(mom_lag_diff)
            #print(mom_obs_diff)

            gamma_tot = gamma_tot .+ mom_obs_diff*mom_lag_diff'
        end
    end

    gamma = (1/(T*H)).*gamma_tot

    return gamma
end


function Bootstrap()
    rho_hat = zeros(100)
    sigma_hat = zeros(100)
        for i = 1:100
            Random.seed!(i)
            P,R = initialize()
            xy= simulate_0(P,R.e_0)
            R.x_t = xy[:,1]
            x_tm1 = zeros(200)
            for t=1:P.T
                if t > 1
                    x_tm1[t] = R.x_t[t-1]
                else
                    x_tm1[t] = 0
                end
            end
            R.M_t = [mean(R.x_t), sum(((R.x_t .- mean(R.x_t)).*((R.x_t .- mean(R.x_t)))))/(P.T*P.H_0), sum(((R.x_t .- mean(R.x_t)).*((x_tm1 .- mean(R.x_t)))))/(P.T*P.H_0)]
            sigma_hat[i],rho_hat[i]  = procedure_smm(P,R,4)
            #println(rho_hat[i])
        end
        return rho_hat, sigma_hat
end
rho_hat,sigma_hat = Bootstrap()
