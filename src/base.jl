export linearLagrangePolynomials, compute1stDeriv4linearLagrangePolynomials, constantFunctions

export constructBasisFunctionMatrixLinearLagrange, constructBasisFunctionMatrixConstantFuncs

export GaussLegendreQuadRule

import GaussQuadrature.legendre


function constructBasisFunctionMatrixLinearLagrange(;interval::AbstractArray=[-1,1], evalPts::AbstractArray=[-1,1])
    N0, N1 = linearLagrangePolynomials(interval=interval)
    dN0, dN1 = compute1stDeriv4linearLagrangePolynomials(interval=interval)

    N_matrix = [(N0.(evalPts))'; (N1.(evalPts))']
    dN_matrix = [(dN0.(evalPts))'; (dN1.(evalPts))']

    return (sparse(N_matrix), sparse(dN_matrix))
end


function constructBasisFunctionMatrixConstantFuncs(;evalPts::AbstractArray=[-1,1])
    R = constantFunctions()
    
    R_matrix = Matrix((R.(evalPts))');

    return sparse(R_matrix)
end


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

