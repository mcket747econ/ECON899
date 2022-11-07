using Parameters, Plots, Printf, DelimitedFiles, DataFrames, Setfield
cd("C:/Users/Rafeh/Documents/GitHub/ECON899/PS6")
include("ps6_model.jl");

## No Shocks
P,R = Initialize()
stationary_price(P,R,0.0)
##Solve for market
R0 = stationary_M(P,R)


## Tauchen Shocks alpha = 1.0
P,R = Initialize()
stationary_price(P,R,1.0)
##Solve for market
R1 = stationary_M(P,R)

## Tauchen Shocks alpha = 2.0
P,R = Initialize()
stationary_price(P,R,2.0)
##Solve for market
R2 = stationary_M(P,R)


## Plot exit share
plot([R0.pf_entry_x R1.pf_entry_x R2.pf_entry_x],
             label = ["No Shocks" "Tauchen Shocks (1.0)" "Tauchen Shocks (2.0)"],
             title = "Decision Rules of Exit")
savefig("plot_regular.png")

## 
data_summary = DataFrame(hcat(["Price","Mass Entering", "Mass Stay", "Mass Exit", "L stayers", "Labor of enter","Agg Labor", "Frac Labor"], 
table(P,R0), table(P,R1), table(P,R2)),["Vals","tv0","tv1","tv2"])

##Save
writedlm("shocks.csv",eachrow(data_summary),",")



## Make cost of firing higher 
P,R = Initialize()
P2 = reconstruct(P,cf=15.0)

## Solve with no shocks
stationary_price(P2,R,0.0)
R0 = stationary_M(P2,R)


## Solve with alpha=1
P,R = Initialize()
stationary_price(P2,R,1.0)
R1 = stationary_M(P2,R)

## Solve with alpha=2
P,R = Initialize()
stationary_price(P2,R,2.0)
R2 = stationary_M(P2,R)


data_summary2 = DataFrame(hcat(["Price","Mass Entering", "Mass Stay", "Mass Exit", "L stayers", "Labor of enter","Agg Labor", "Frac Labor"], 
table(P,R0), table(P,R1), table(P,R2)),["Vals","tv0","tv1","tv2"])

## Export
writedlm("higher_cf.csv",eachrow(data_summary2),",")


## Plot
plot([R0.pf_entry_x R1.pf_entry_x R2.pf_entry_x],
             label = ["No Shocks" "Tauchen Shocks (1.0)" "Tauchen Shocks (2.0)"],
             title = "Decision Rules of Exit")
savefig("plot_higher_cf.png")
