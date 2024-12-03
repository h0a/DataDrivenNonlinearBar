export linearLagrangePolynomials, compute1stDeriv4linearLagrangePolynomials, constantFunctions


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

