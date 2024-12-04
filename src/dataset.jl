export generateDataHookLaw, generateDataSigmoid

import Noise.add_gauss


function generateDataHookLaw(;Young_modulus::Float64=1.0, numDataPts::Int=100, strainRange::AbstractArray=[0,1], add_noise::Bool=false, noise_magnitude::Float64=0.1)

    # input:
    # noise_magnitude: percentage of the maximum stress value as added noise
    
    e0, e1 = strainRange
    strain_data = collect(LinRange(e0,e1,numDataPts))
    stress_data = Young_modulus .* strain_data

    if add_noise
        stress_data = add_gauss(stress_data, noise_magnitude*maximum(stress_data));
    end

    return [strain_data stress_data]
end





function generateDataSigmoid(;max_stress::Float64=1.0, numDataPts::Int=100, strainRange::AbstractArray=[0,1], add_noise::Bool=false, noise_magnitude::Float64=0.1)

    # input:
    # noise_magnitude: percentage of the maximum stress value as added noise
    
    e0, e1 = strainRange
    strain_data = collect(LinRange(e0,e1,numDataPts))
    stress_data = max_stress .* ( 2 ./ (1 .+ exp.(-strain_data .* 1e3)) .- 1 )

    if add_noise
        stress_data = add_gauss(stress_data, noise_magnitude*maximum(stress_data));
    end

    return [strain_data stress_data]
end