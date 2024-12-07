export assembleLinearSystemMatrix, assembleRhsLinearBar
export assembleBalanceResidual, assembleLinearizedSystemMatrix


function assembleRhsLinearBar(;data_star::AbstractArray, node_vector::AbstractArray, num_ele::Int=2, numQuadPts::Int=2, ndofs::AbstractArray=ones(Int,5), costFunc_constant::Float64=1.0, bar_distF::Float64=1.0, cross_section_area::Float64=1.0)

    # quad points in default interval [-1,1]
    quad_pts, quad_weights = GaussLegendreQuadRule(numQuadPts=numQuadPts);
    
    # basis function matrix evaluated in master element [-1,1]
    N_matrix, dN_matrix = constructBasisFunctionMatrixLinearLagrange(evalPts=quad_pts)
    R_matrix = constructBasisFunctionMatrixConstantFuncs(evalPts=quad_pts)

    # allocation blocks of rhs
    ndof_u, ndof_e, ndof_s, ndof_mu, ndof_lambda = ndofs;
    rhs_b1 = spzeros(ndof_u);
    rhs_b5 = spzeros(ndof_lambda)

    rhs_b2 = spzeros(ndof_e);
    rhs_b3 = spzeros(ndof_s)
    rhs_b4 = spzeros(ndof_mu)

    # assembly routine
    for cc_ele = 1:num_ele      # loop over elements    
        active_dofs_u = active_dofs_lambda = collect(cc_ele:cc_ele+1)
        active_dofs_e = active_dofs_s = active_dofs_mu = cc_ele
    
        # jacobian for the integration
        xi0, xi1 = node_vector[cc_ele:cc_ele+1];
        J4int = (xi1 - xi0) / 2;   
       
        # integrated blocks of the rhs
        rhs_b2[active_dofs_e] += ( R_matrix * (cross_section_area .* quad_weights .* J4int .* costFunc_constant .* (R_matrix' * data_star[active_dofs_e,1])) )[1]
        
        rhs_b3[active_dofs_s] += ( R_matrix * (cross_section_area .* quad_weights .* J4int ./ costFunc_constant .* (R_matrix' * data_star[active_dofs_s,2])) )[1]
    
        rhs_b5[active_dofs_lambda] += N_matrix * (quad_weights .* J4int .* bar_distF)
    end    
    
    # global rhs
    rhs = [rhs_b1; rhs_b2; rhs_b3; rhs_b4; rhs_b5];

    return rhs
end



function assembleLinearSystemMatrix(;node_vector::AbstractArray, num_ele::Int=2, numQuadPts::Int=2, ndofs::AbstractArray=ones(Int,5), costFunc_constant::Float64=1.0, cross_section_area::Float64=1.0)

    # quad points in default interval [-1,1]
    quad_pts, quad_weights = GaussLegendreQuadRule(numQuadPts=numQuadPts);
    
    # basis function matrix evaluated in master element [-1,1]
    N_matrix, dN_matrix = constructBasisFunctionMatrixLinearLagrange(evalPts=quad_pts)
    R_matrix = constructBasisFunctionMatrixConstantFuncs(evalPts=quad_pts)
    
    # alloccation blocks of the system matrix
    ndof_u, ndof_e, ndof_s, ndof_mu, ndof_lambda = ndofs;
    J11, J12, J13, J14, J15 = spzeros(ndof_u,ndof_u), spzeros(ndof_u,ndof_e), spzeros(ndof_u,ndof_s), spzeros(ndof_u,ndof_mu), spzeros(ndof_u,ndof_lambda);

    J22, J23, J24, J25 = spzeros(ndof_e,ndof_e), spzeros(ndof_e,ndof_s), spzeros(ndof_e,ndof_mu), spzeros(ndof_e,ndof_lambda);

    J33, J34, J35 = spzeros(ndof_s,ndof_s), spzeros(ndof_s,ndof_mu), spzeros(ndof_s,ndof_lambda);

    J44, J45 = spzeros(ndof_mu,ndof_mu), spzeros(ndof_mu,ndof_lambda);

    J55 = spzeros(ndof_lambda,ndof_lambda);

    # assembly routine
    for cc_ele = 1:num_ele      # loop over elements
        active_dofs_u = active_dofs_lambda = collect(cc_ele:cc_ele+1)
        active_dofs_e = active_dofs_s = active_dofs_mu = cc_ele

        # jacobian for the integration
        xi0, xi1 = node_vector[cc_ele:cc_ele+1];
        J4int = (xi1 - xi0) / 2;

        # jacobian for derivative
        J4deriv = (xi1 - xi0) / 2;

        # integrated blocks of the system matrix        
        J14[active_dofs_u,active_dofs_mu] += (dN_matrix ./ J4deriv) * (cross_section_area .* quad_weights .* J4int .* R_matrix')
        
        J22[active_dofs_e,active_dofs_e] += (R_matrix * ( cross_section_area .* quad_weights .* J4int .* costFunc_constant .* R_matrix' ))[1]

        J24[active_dofs_e,active_dofs_mu] += -(R_matrix * ( cross_section_area .* quad_weights .* J4int .* R_matrix' ))[1]

        J33[active_dofs_s,active_dofs_s] += (R_matrix * ( cross_section_area .* quad_weights .* J4int ./ costFunc_constant .* R_matrix' ))[1]

        J35[active_dofs_s,active_dofs_lambda] += (R_matrix * (cross_section_area .* quad_weights .* J4int .* (dN_matrix ./ J4deriv)'))[:]
    end


    # global system matrix
    J = [J11  J12  J13  J14  J15;
         J12' J22  J23  J24  J25;
         J13' J23' J33  J34  J35;
         J14' J24' J34' J44  J45;
         J15' J25' J35' J45' J55];

    return J
end



function assembleBalanceResidual(;x::AbstractArray, data_star::AbstractArray, node_vector::AbstractArray, num_ele::Int=2, numQuadPts::Int=2, ndofs::AbstractArray=ones(Int,5), costFunc_constant::Float64=1.0, bar_distF::Float64=1.0, cross_section_area::Float64=1.0)

    # quad points in default interval [-1,1]
    quad_pts, quad_weights = GaussLegendreQuadRule(numQuadPts=numQuadPts);
    
    # basis function matrix evaluated in master element [-1,1]
    N_matrix, dN_matrix = constructBasisFunctionMatrixLinearLagrange(evalPts=quad_pts)
    R_matrix = constructBasisFunctionMatrixConstantFuncs(evalPts=quad_pts)

    # extract variable fields from solution vector of the previous iteration
    ndof_u, ndof_e, ndof_s, ndof_mu, ndof_lambda = ndofs;

    uhat = x[1:ndof_u]
    ebar = x[ndof_u+1:ndof_u+ndof_e]
    sbar = x[ndof_u+ndof_e+1:ndof_u+ndof_e+ndof_s]
    mubar = x[ndof_u+ndof_e+ndof_s+1:ndof_u+ndof_e+ndof_s+ndof_mu]
    lambdahat = x[ndof_u+ndof_e+ndof_s+ndof_mu+1:end]

    # prepare the difference between material data
    e_diff = ebar - data_star[:,1]
    s_diff = sbar - data_star[:,2]

    # alloccation blocks of rhs
    rhs_b1 = spzeros(ndof_u)
    rhs_b5 = spzeros(ndof_lambda)

    rhs_b2 = spzeros(ndof_e)
    rhs_b3 = spzeros(ndof_s)
    rhs_b4 = spzeros(ndof_mu)

    # assembly routine
    for cc_ele = 1:num_ele      # loop over elements    
        active_dofs_u = active_dofs_lambda = collect(cc_ele:cc_ele+1)
        active_dofs_e = active_dofs_s = active_dofs_mu = cc_ele
    
        # jacobian for the integration
        xi0, xi1 = node_vector[cc_ele:cc_ele+1];
        J4int = (xi1 - xi0) / 2;
    
        # jacobian for derivative
        J4deriv = (xi1 - xi0) / 2;
    
        # interpolate discrete solution from the previous iteration
        duh = (dN_matrix ./ J4deriv)' * uhat[active_dofs_u]
    
        eh = R_matrix' * ebar[active_dofs_e]
        sh = R_matrix' * sbar[active_dofs_s]
        muh = R_matrix' * mubar[active_dofs_mu]
    
        dlambdah = (dN_matrix ./ J4deriv)' * lambdahat[active_dofs_lambda]
    
        # integrated blocks of the rhs
        rhs_b1[active_dofs_u] += - (dN_matrix ./ J4deriv) * (cross_section_area .* quad_weights .* J4int .* (sh .* dlambdah)) - 
                                   (dN_matrix ./ J4deriv) * (cross_section_area .* quad_weights .* J4int .* ((1 .+ duh) .* muh))
        
        rhs_b2[active_dofs_e] += ( R_matrix * (cross_section_area .* quad_weights .* J4int .* muh) - 
                                   R_matrix * (cross_section_area .* quad_weights .* J4int .* costFunc_constant .* (R_matrix' * e_diff[active_dofs_e])) )[1]
        
        rhs_b3[active_dofs_s] += ( -R_matrix * (cross_section_area .* quad_weights .* J4int .* (1 .+ duh) .* dlambdah) - 
                                    R_matrix * (cross_section_area .* quad_weights .* J4int ./ costFunc_constant .* (R_matrix' * s_diff[active_dofs_s])))[1]
    
        rhs_b4[active_dofs_mu] += ( -R_matrix * (cross_section_area .* quad_weights .* J4int .* (duh + 0.5 .* duh .* duh - eh)))[1]
    
        rhs_b5[active_dofs_lambda] += N_matrix * (quad_weights .* J4int .* bar_distF) - 
                                     (dN_matrix ./ J4deriv) * (cross_section_area .* quad_weights .* J4int .* ((1 .+ duh) .* sh))
    end    
    
    # global rhs
    rhs = [rhs_b1; rhs_b2; rhs_b3; rhs_b4; rhs_b5];

    return rhs
end



function assembleLinearizedSystemMatrix(;x::AbstractArray, node_vector::AbstractArray, num_ele::Int=2, numQuadPts::Int=2, ndofs::AbstractArray=ones(Int,5), costFunc_constant::Float64=1.0, cross_section_area::Float64=1.0)

    # quad points in default interval [-1,1]
    quad_pts, quad_weights = GaussLegendreQuadRule(numQuadPts=numQuadPts);
    
    # basis function matrix evaluated in master element [-1,1]
    N_matrix, dN_matrix = constructBasisFunctionMatrixLinearLagrange(evalPts=quad_pts)
    R_matrix = constructBasisFunctionMatrixConstantFuncs(evalPts=quad_pts)


    # extract variable fields from solution vector of the previous iteration
    ndof_u, ndof_e, ndof_s, ndof_mu, ndof_lambda = ndofs;

    uhat = x[1:ndof_u]
    ebar = x[ndof_u+1:ndof_u+ndof_e]
    sbar = x[ndof_u+ndof_e+1:ndof_u+ndof_e+ndof_s]
    mubar = x[ndof_u+ndof_e+ndof_s+1:ndof_u+ndof_e+ndof_s+ndof_mu]
    lambdahat = x[ndof_u+ndof_e+ndof_s+ndof_mu+1:end]


    # alloccation blocks of the system matrix
    J11, J12, J13, J14, J15 = spzeros(ndof_u,ndof_u), spzeros(ndof_u,ndof_e), spzeros(ndof_u,ndof_s), spzeros(ndof_u,ndof_mu), spzeros(ndof_u,ndof_lambda);

    J22, J23, J24, J25 = spzeros(ndof_e,ndof_e), spzeros(ndof_e,ndof_s), spzeros(ndof_e,ndof_mu), spzeros(ndof_e,ndof_lambda);

    J33, J34, J35 = spzeros(ndof_s,ndof_s), spzeros(ndof_s,ndof_mu), spzeros(ndof_s,ndof_lambda);

    J44, J45 = spzeros(ndof_mu,ndof_mu), spzeros(ndof_mu,ndof_lambda);

    J55 = spzeros(ndof_lambda,ndof_lambda);


    # assembly routine
    for cc_ele = 1:num_ele      # loop over elements
        active_dofs_u = active_dofs_lambda = collect(cc_ele:cc_ele+1)
        active_dofs_e = active_dofs_s = active_dofs_mu = cc_ele

        # jacobian for the integration
        xi0, xi1 = node_vector[cc_ele:cc_ele+1];
        J4int = (xi1 - xi0) / 2;

        # jacobian for derivative
        J4deriv = (xi1 - xi0) / 2;

        # interpolate discrete solution from the previous iteration
        duh = (dN_matrix ./ J4deriv)' * uhat[active_dofs_u]

        eh = R_matrix' * ebar[active_dofs_e]
        sh = R_matrix' * sbar[active_dofs_s]
        muh = R_matrix' * mubar[active_dofs_mu]

        dlambdah = (dN_matrix ./ J4deriv)' * lambdahat[active_dofs_lambda]

        # integrated blocks of the system matrix
        J11[active_dofs_u,active_dofs_u] += (dN_matrix ./ J4deriv) * (cross_section_area .* quad_weights .* J4int .* muh .* (dN_matrix ./ J4deriv)')
        
        J13[active_dofs_u,active_dofs_s] += (dN_matrix ./ J4deriv) * (cross_section_area .* quad_weights .* J4int .* dlambdah .* R_matrix')
        
        J14[active_dofs_u,active_dofs_mu] += (dN_matrix ./ J4deriv) * (cross_section_area .* quad_weights .* J4int .* ((1 .+ duh) .* R_matrix'))
        
        J15[active_dofs_u,active_dofs_lambda] += (dN_matrix ./ J4deriv) * (cross_section_area .* quad_weights .* J4int .* sh .* (dN_matrix ./ J4deriv)')
        
        J22[active_dofs_e,active_dofs_e] += ( R_matrix * (cross_section_area .* quad_weights .* J4int .* costFunc_constant .* R_matrix') )[1]

        J24[active_dofs_e,active_dofs_mu] += -( R_matrix * (cross_section_area .* quad_weights .* J4int .* R_matrix') )[1]

        J33[active_dofs_s,active_dofs_s] += ( R_matrix * (cross_section_area .* quad_weights .* J4int ./ costFunc_constant .* R_matrix') )[1]

        J35[active_dofs_s,active_dofs_lambda] += ( R_matrix * (cross_section_area .* quad_weights .* J4int .* (1 .+ duh) .* (dN_matrix ./ J4deriv)') )[:]
    end


    # global system matrix
    J = [J11  J12  J13  J14  J15;
         J12' J22  J23  J24  J25;
         J13' J23' J33  J34  J35;
         J14' J24' J34' J44  J45;
         J15' J25' J35' J45' J55];

    return J
end