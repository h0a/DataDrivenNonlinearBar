export directSolverLinearBar


function directSolverLinearBar(;node_vector::AbstractArray, data_set::AbstractArray, costFunc_ele::Function, random_init_data::Bool=true, max_iter::Int=2000, tol::Float64=1e-10, num_ele::Int=2, numQuadPts::Int=2, costFunc_constant::Float64=1.0, bar_distF::Float64=1.0)

    ## initialize e_star and s_star
    numDataPts = size(data_set,1)

    if random_init_data
        init_data_id = rand(1:numDataPts,num_ele)
        init_data = data_set[init_data_id,:]
    end

    data_star = deepcopy(init_data)

    ## direct solver [Kirchdoerfer - 2016 - Data-driven computational mechanics]
    # ndofs
    num_node = length(node_vector)
    ndof_u = ndof_lambda = num_node
    ndof_e = ndof_s = ndof_mu = num_ele

    ndof_tot = ndof_u + ndof_e + ndof_s + ndof_mu + ndof_lambda
    ndofs = [ndof_u, ndof_e, ndof_s, ndof_mu, ndof_lambda]

    # assembly system matrix
    A = assembleLinearSystemMatrix(node_vector=node_vector, num_ele=num_ele, ndofs=ndofs, costFunc_constant=costFunc_constant)

    # boundary conditions: fixed-free
    constrained_dofs = [1 (ndof_u+ndof_e+ndof_s+ndof_mu+1)]
    ids = collect(1:size(A,1))
    deleteat!(ids, constrained_dofs)
    A = A[ids,ids]

    # check condition number
    κ = cond(Matrix(A))
    @show κ

    # allocation solution vector
    x = spzeros(ndof_tot-length(constrained_dofs))
    costFunc_global = []

    # iterative direct solver
    iter = 0

    while iter <= max_iter
        # assembly rhs
        rhs = assembleRhsLinearBar(data_star=data_star, node_vector=node_vector, num_ele=num_ele, ndofs=ndofs, costFunc_constant=costFunc_constant, bar_distF=bar_distF)
        rhs = rhs[ids]

        # solving
        x = qr(Matrix(A)) \ rhs

        # check residual (Ax-rhs)
        r = norm(rhs - A * x)
        @show r

        # collect computed ebar and sbar
        ebar = x[ndof_u:ndof_u+ndof_e-1]
        sbar = x[ndof_u+ndof_e:ndof_u+ndof_e+ndof_s-1]

        ## local state assignment
        data_star_new = assignLocalState(data_set=data_set, local_state=[ebar sbar], costFunc_ele=costFunc_ele)

        # evaluate global cost function (discrete)
        push!(costFunc_global, integrateCostfunction(costFunc_ele=costFunc_ele, local_state=[ebar sbar], data_star=data_star, node_vector=node_vector, num_ele=num_ele, numQuadPts=numQuadPts))

        ## test convergence
        data_diff = data_star - data_star_new
        err = [ norm(data_diff[i,:]) for i in 1:num_ele ]
        max_err = maximum(err)

        iter += 1

        if max_err <= tol
            @show iter
            break
        end

        # overwrite local state
        data_star = deepcopy(data_star_new)        
    end

    # collect solution fields
    uhat, ebar, sbar = [0; x[1:ndof_u-1]], x[ndof_u:ndof_u+ndof_e-1], x[ndof_u+ndof_e:ndof_u+ndof_e+ndof_s-1]

    return uhat, ebar, sbar, costFunc_global
end



function assignLocalState(;data_set::AbstractArray, local_state::AbstractArray, costFunc_ele::Function)
    num_local_state = size(local_state,1)
    numDataPts = size(data_set,1)

    # allocation
    data_star = zeros(num_local_state,2)
    
    # find the closest data point to the local state
    for ii = 1:num_local_state
        d = zeros(numDataPts);
    
        for i = 1:numDataPts
            d[i] = sqrt( costFunc_ele(local_state[ii,1] - data_set[i,1], local_state[ii,2] - data_set[i,2]) );
        end
        min_id = findmin(d)[2];
        data_star[ii,:] = data_set[min_id,:];
    end

    return data_star
end



function integrateCostfunction(;costFunc_ele::Function, local_state::AbstractArray, data_star::AbstractArray, node_vector::AbstractArray, num_ele::Int=2, numQuadPts::Int=2)

    # quad points in default interval [-1,1]
    quad_pts, quad_weights = GaussLegendreQuadRule(numQuadPts=numQuadPts);

    # integration
    costFunc_global = 0

    for e in 1:num_ele      # loop over element
        # jacobian for the integration
        xi0, xi1 = node_vector[e:e+1]
        J4int = (xi1 - xi0) / 2

        costFunc_global += costFunc_ele(local_state[e,1] - data_star[e,1], local_state[e,2] - data_star[e,2]) * sum(quad_weights) * J4int
    end

    return costFunc_global
end
