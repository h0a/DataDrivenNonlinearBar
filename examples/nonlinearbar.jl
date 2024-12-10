
using LinearAlgebra, SparseArrays, StaticArrays, Statistics, Plots
using DataDrivenNonlinearBar


#### data-driven for nonlinear fixed-free bar with constant uniform distributed load and linear material law


# inputs
bar_len = 2.0;      # [m]   - initial length of the bar
bar_area = 1.5;     # [m^2] - cross-sectional area of the bar
bar_distF = 1.8;    # [N]   - constant uniform distributed load
bar_E = 1e4;        # [Pa]  - assumed Young_modulus
num_ele = 10;       # [-]   - number of elements
numDataPts = 200;   # [-]   - number of data points


# generate data: linear function
strain_limit = 2*bar_distF * bar_len / (bar_E * bar_area);
SE = generateDataHookLaw(Young_modulus=bar_E, numDataPts=numDataPts, strainRange=[-strain_limit,strain_limit]);


# elementwise cost function
costFunc_constant =  mean(SE[:,2] ./ SE[:,1]);
costFunc_ele = (e,s) -> 0.5 * (costFunc_constant * e^2 + 1/costFunc_constant * s^2);


# node vector
num_node = num_ele + 1;
node_vector = collect(LinRange(0.0,bar_len,num_node));


# solving
uhat, ebar, sbar, costFunc_global = directSolverNonLinearBar(node_vector=node_vector, data_set=SE, costFunc_ele=costFunc_ele, num_ele=num_ele, costFunc_constant=costFunc_constant, bar_distF=bar_distF, cross_section_area=bar_area);


# plots
plot(SE[:,1], SE[:,2])
scatter!(ebar,sbar)


plot(1:length(costFunc_global), costFunc_global, xscale=:log10, yscale=:log10)

