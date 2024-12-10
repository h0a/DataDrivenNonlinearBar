export NewtonRaphsonStep


function NewtonRaphsonStep(;previous_sol::AbstractArray, data_star::AbstractArray, node_vector::AbstractArray, num_ele::Int=2, numQuadPts::Int=2, ndofs::AbstractArray=ones(Int,5), costFunc_constant::Float64=1.0, bar_distF::Float64=1.0, cross_section_area::Float64=1.0, constrained_dofs::AbstractArray=[1])

    ndof_u, ndof_e, ndof_s, ndof_mu, ndof_lambda = ndofs;

    # assembly
   rhs = assembleBalanceResidual(x=previous_sol, data_star=data_star, node_vector=node_vector, num_ele=num_ele, numQuadPts=numQuadPts, ndofs=ndofs, costFunc_constant=costFunc_constant, bar_distF=bar_distF, cross_section_area=cross_section_area);

    J = assembleLinearizedSystemMatrix(x=previous_sol, node_vector=node_vector, num_ele=num_ele, numQuadPts=numQuadPts, ndofs=ndofs, costFunc_constant=costFunc_constant, cross_section_area=cross_section_area);

    # enforcing boundary conditions    
    ids = collect(1:size(J,1))
    deleteat!(ids, constrained_dofs)

    J = J[ids,ids];
    rhs = rhs[ids];

    # solving
    Delta_x = qr(Matrix(J)) \ rhs
    # Delta_x = Matrix(J) \ rhs

    # check residual and condition number of J
    r = norm(rhs - J * Delta_x)
    κ = cond(Matrix(J))    
    
    if r > 1e-10
        # @show r
        print("Warning: Solution with residual > 1e-10")
    end
    if κ > 1e20
        # @show κ
        print("Warning: Condition number of the system matrix > 1e20")
    end
    
    return Delta_x
end