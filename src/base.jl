export linearLagrangePolynomials, compute1stDeriv4linearLagrangePolynomials, constantFunctions

export GaussLegendreQuadRule

import GaussQuadrature.legendre

# linear Lagrange polynomial on an interval [x0,x1]
function linearLagrangePolynomials(;interval::AbstractArray=[-1,1])
    x0, x1 = interval
    L0 = x -> (x-x1) / (x0-x1)
    L1 = x -> (x-x0) / (x1-x0)
    
    return L0,L1
end


function compute1stDeriv4linearLagrangePolynomials(;interval::AbstractArray=[-1,1])
    x0, x1 = interval
    dL0 = x -> 1 / (x0-x1)
    dL1 = x -> 1 / (x1-x0)
    
    return dL0,dL1
end


function constantFunctions()    
    L = x -> 1.0
    
    return L
end


# Gauss-Legendre quadrature rule

function GaussLegendreQuadRule(;interval::AbstractArray=[-1,1], numQuadPts::Int=2)
    # Gauss-Legendre quadrature rule on [-1,1]
    x, w = legendre(numQuadPts)

    # map to an arbitrary interval
    x0,x1 = interval
    m, n = (x1-x0)/2, (x1+x0)/2

    return x .* m .+ n, w .* m
end

