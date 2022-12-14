using StatFiles, DataFrames, Optim, LinearAlgebra, LatexPrint, ForwardDiff, CSV

cd("/Users/jacobbills/Desktop/Economics/Econ 899/JF 1/")
include("ps1b_model.jl") #source of functions

data = DataFrame(load("Mortgage_performance_data.dta"))

y = data[!,20]
y = Array(y)
ind = ["i_large_loan", "i_medium_loan" , "rate_spread" , "i_refinance", "age_r", "cltv", "dti", "cu" , "first_mort_r", "score_0", "score_1", "i_FHA", "i_open_year2", "i_open_year3", "i_open_year4", "i_open_year5"]
x = data[!,ind]
insertcols!(x, 17, :const => 1) 
x = Matrix(x)
β::Array{Float64,1} = zeros(size(x,2)) #I counted sixteen variables in the problem set, plus a constant
β[17] = -1 #if I understand the problem right


#Question 1

ll1, score1, H1 = evaluate_b(β,y,x)

println("The Log-Likelihood is ", ll1)
println("The score of the LL is: ", score1)
println("The Hessian of the LL is: ", H1)


lap(score1) #LaTeX printout of score
tabular(round.(H1,digits=2)) #LaTeX printout of Hessian

##Question 2 (not sure how to approach this)


g, h = q2(β,y,x) #function defined to do question 2 

println(norm(g-score1)) #how different are the vectors? Not very! 
println(norm(h-H1))




#Question 3 
@time b3, ll3= newton(β, y, x)

println("The coeffcients are: ", b3)
println("The log-likelihood is ", ll3)

#Question 4

fun(b) = -LL(b,y,x) #define a new function for the purposes of minimizing only β
@time opta = optimize(fun, β, BFGS()) #optimize the function using BFGS 
b4a = Optim.minimizer(opta) #get coefficients 
println("The BFGS coefficients are: ", b4a)

@time optb = optimize(fun, β, NelderMead(), Optim.Options(iterations=5000)) #Simplex method, up-ed the number of iterations so it wouldn't fail 
b4b = Optim.minimizer(optb) #get coefficients 
println("The Simplex coefficients are: ", b4b)

println(norm(b4a-b3))

println(norm(b4b-b3))

#Print Results
CSV.write("q1_score.csv", DataFrame([score1],[:col1])) 
CSV.write("q1_hessian.csv",DataFrame(H1, :auto))

CSV.write("q2_first.csv", DataFrame([g],[:col1])) 
CSV.write("q2_second.csv", DataFrame(h, :auto))

CSV.write("q3_coef.csv", DataFrame([b3],[:col1]))

CSV.write("q4_BFGS.csv", DataFrame([b4a],[:col1]))
CSV.write("q4_Simplex.csv", DataFrame([b4b],[:col1]))