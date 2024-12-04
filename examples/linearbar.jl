
using LinearAlgebra, SparseArrays, StaticArrays, Statistics
using DataDrivenNonlinearBar


# generate data: either linear or sigmoid function
numDataPts = 200;
SE = generateDataHookLaw(Young_modulus=1e4, numDataPts=numDataPts, strainRange=[-2e-4,2e-4]);

# elementwise cost function
costFunc_constant =  mean(SE[:,2] ./ SE[:,1]);
costFunc_ele = (e,s) -> 0.5 * (costFunc_constant * e^2 + 1/costFunc_constant * s^2);


# parameters
bar_len = 1.0;      # [m] - initial length of the bar
bar_area = 2e-3;    # [m^2] - cross-sectional area of the bar
bar_distF = 1.0;    # [N] - constant uniform distributed load

# node vector
num_ele = 10;
num_node = num_ele + 1;
node_vector = collect(LinRange(0.0,bar_len,num_node));

# solving
uhat, ebar, sbar, costFunc_global = directSolverLinearBar(node_vector=node_vector, data_set=SE, costFunc_ele=costFunc_ele, num_ele=num_ele, costFunc_constant=costFunc_constant, bar_distF=bar_distF);

# plots
using Plots

plot(SE[:,1], SE[:,2])
scatter!(ebar,sbar)


plot(1:length(costFunc_global), costFunc_global, xscale=:log10, yscale=:log10)
