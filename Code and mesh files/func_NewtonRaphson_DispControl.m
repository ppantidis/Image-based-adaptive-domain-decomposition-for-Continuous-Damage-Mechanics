function [Domain_Vector, Res_F_F_norm, Res_u_norm, iteration,Reactions_y,Khh_inv,lr_vec_inclist,Beta] ...
    = func_NewtonRaphson_DispControl(fixnodes_applied,increment,Domain_Vector, K_pos_mat,ndomains,IsProj,lr,Khh_inv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ==================== NEWTON RAPHSON ANALYSIS METHOD =====================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include global variables
func_include_flags;
schurcomplement = 1;

% Find appropriate penalty parameter value
Beta = 0;
for i = 1:1
    Beta = max(max(max(Domain_Vector(i).K)),Beta);
end
Beta = Beta * 1E4;
% Beta = max(max(Domain_Vector(1).K)) * 1e4;

% tolerance below which solution should stop
tol = 1*10^-5;

% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOLVER ROUTINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
if IsProj == 0
    
    for i = 1:ndomains
        % Evaluate stiffness matrix [K] and history variable matrix at the start of the increment
        Domain_Vector(i).history_var_mat_previousinc = Domain_Vector(i).history_var_mat_stored;
        if SolverID == 3 && Domain_Vector(i).DomainDamage
            Domain_Vector(i).damage_mat_inl_prev_inc = Domain_Vector(i).damage_mat_inl_prev_inc_stored;
        end
        Domain_Vector(i) = func_domain_global_stiffness(Domain_Vector(i),IsProj);
    end

    % Initiate the iteration counter within each loading step
    iteration = 1;
    [K,DRdu,dofs,~,~,nnodes] = func_assemble(Domain_Vector,ndomains, K_pos_mat,Beta);

    % Partition K and DRdu
    if TangentID == 1
        [~, ~, ~, K_F, ~, ~, listofnodes_ebc, listofnodes_nbc] = func_partitionK(Domain_Vector(1).ndof,nnodes,fixnodes_applied, K);
    elseif TangentID == 2
        [~, ~, ~, DRdu_F, ~, ~, listofnodes_ebc, listofnodes_nbc] = func_partitionK(Domain_Vector(1).ndof,nnodes,fixnodes_applied, DRdu);
    else
        disp("Check your choice of solver!")
    end

    % Calculate the fixed displacement increment and partition dofs
    [~, dofs_E, dofs_F] = func_partitiond_fixnodes(Domain_Vector(1).ndof,nnodes,fixnodes_applied,dofs);

    % Assemble the global dofs
    dofs(listofnodes_ebc) = dofs_E;
    dofs(listofnodes_nbc) = dofs_F;

    position = 0;
    for i = 1:ndomains
        Domain_Vector(i).dofs = dofs(position+1:position+length(Domain_Vector(i).dofs));
        position              = position + length(Domain_Vector(i).dofs);
    end

    for i = 1:ndomains
        % Evaluate the updated stiffness matrix [K] and history variable matrix
%         tic
        Domain_Vector(i) = func_domain_global_stiffness(Domain_Vector(i),IsProj);
%         toc
    end

    [~, ~, dofs, Res_F, f_external, nnodes] = func_assemble(Domain_Vector,ndomains,K_pos_mat,Beta);

    % Calculate and partition the residual force vector (based on internal stresses) after delta_dofs_E has been applied
    [~, ~, Res_F_F] = func_partitionf(Domain_Vector(1).ndof,nnodes,fixnodes_applied,Res_F);
    [~, ~, f_ext_F] = func_partitionf(Domain_Vector(1).ndof,nnodes,fixnodes_applied,f_external);

    % Compute and store the norm of the residual at the free boundary
    Res_F_F_norm_first      = norm(Res_F_F,2);
    Res_F_F_norm(iteration) = Res_F_F_norm_first;

    % The while loop runs until either the error is smaller than a specific
    % threshold (convergence is achieved) or the number of iterations within
    % a loading step is exceeded (computational cost is high)
    while true
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate the displacements increment at the natural boundary
        if schurcomplement == 0 || ndomains == 1
%             disp("First block:")
%             tic
            if TangentID == 1
                if  iteration > 2
                    if Res_u_norm(iteration-1) > mean(Res_u_norm(max(1,iteration-2):iteration-2))
                        lr = min(lr * 10, 1e-3); % lr * 10
                    else
                        lr = lr / 10;
                    end
                end
                delta_dofs_F = - (K_F + lr*diag(diag(K_F))) \ Res_F_F;     % Analytical calculation (using K)
                % delta_dofs_F = - sparse(K_F + lr*diag(diag(K_F))) \ Res_F_F;     % Analytical calculation (using K)                
            elseif TangentID == 2
                delta_dofs_F = - DRdu_F \ Res_F_F;  % Numerical approximation (using DRdu)
            else
                disp("Check your choice of solver!")
            end
%             toc
        else
            if(iteration == 1 && newDecomp == 1)
%                 disp("Second block - a")
%                 tic
                if TangentID == 1
                    [Khh_inv,CA_inverse,CA_inverseB,lenA,lenD] = func_getSchurBlocks((K_F+lr*diag(diag(K_F))), Domain_Vector,listofnodes_ebc);
                    delta_dofs_F = func_deltaU(Khh_inv,CA_inverse,CA_inverseB,(K_F+lr*diag(diag(K_F))),Res_F_F,lenA,lenD);
                    newDecomp    = 0;
                else
                    [Khh_inv,CA_inverse,CA_inverseB,lenA,lenD] = func_getSchurBlocks(DRdu_F, Domain_Vector,listofnodes_ebc);
                    delta_dofs_F = func_deltaU(Khh_inv,CA_inverse,CA_inverseB,DRdu_F,Res_F_F,lenA,lenD);
                    newDecomp    = 0;
                end
%                 toc
            else
%                 disp("Second block - b")        
%                 tic
                if TangentID == 1
                    if iteration > 2
                        if Res_u_norm(iteration-1) > mean(Res_u_norm(max(1,iteration-2):iteration-2))
                            lr = min(lr * 10, 1e-3); % lr * 10
                        else
                            lr = lr / 10;
                        end
                    end
                    delta_dofs_F = func_deltaU(Khh_inv,CA_inverse,CA_inverseB,(K_F+lr*diag(diag(K_F))),Res_F_F,lenA,lenD);
                elseif TangentID == 2
                    delta_dofs_F = func_deltaU(Khh_inv,CA_inverse,CA_inverseB,DRdu_F,Res_F_F,lenA,lenD);
                end
%                 toc
            end
        end
        
        lr_vec_inclist(1,iteration) = lr;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update the displacements at the natural boundary
%         disp("Third block")
%         tic
        dofs_F = dofs_F + delta_dofs_F;

        % Assemble the global dofs vector
        dofs(listofnodes_ebc) = dofs_E;
        dofs(listofnodes_nbc) = dofs_F;

        position = 0;
        for i = 1:ndomains
            Domain_Vector(i).dofs = dofs(position+1:position+length(Domain_Vector(i).dofs));
            position              = position + length(Domain_Vector(i).dofs);
        end
        
        for i = 1:ndomains
            % Evaluate stiffness matrix [K] and history variable matrix at the start of the increment
            if Domain_Vector(i).DomainDamage == 1
%                 disp("Hey2")
%                 tic
                Domain_Vector(i) = func_domain_global_stiffness(Domain_Vector(i),IsProj);
%                 toc
            else
                % If DamageDomain is zero, then K will not change, so use the K you already have
                Domain_Vector(i).Res_F = Domain_Vector(i).K*Domain_Vector(i).dofs; 
            end
        end
        
        [K, DRdu, dofs, Res_F, f_external, nnodes] = func_assemble(Domain_Vector,ndomains, K_pos_mat,Beta);

        % Partition the stiffness matrix [K]/[DRdu] and residual force vector {Res_F}
        if TangentID == 1
            [K_E, K_EF, ~, K_F, ~, ~, ~, ~] = func_partitionK(Domain_Vector(1).ndof,nnodes,fixnodes_applied,K);
        elseif TangentID == 2
            [DRdu_E, DRdu_EF, ~, DRdu_F, ~, ~, ~, ~] = func_partitionK(Domain_Vector(1).ndof,nnodes,fixnodes_applied,DRdu);
        else
            disp("Check your choice of solver!")
        end

        [~, Res_F_u, Res_F_F] = func_partitionf(Domain_Vector(1).ndof,nnodes,fixnodes_applied,Res_F);

        % Calculate the norm of the residuals based on a) internal stresses (Res_F_F_norm) and b) displacements (delta_dofs_F):
        Res_F_F_norm(iteration+1) = norm(Res_F_F,2);
        Res_u_norm(iteration)     = norm(delta_dofs_F,2);

        % Keep track of the residual at the first iteration
        if iteration == 1
            Res_u_norm_first = Res_u_norm(iteration);
        end    
        
%         toc

        % Convergence check: 
        if (Res_u_norm(iteration) < tol) || (iteration == max_accept_iter)        
%         if (Res_u_norm(iteration) / Res_u_norm_first  < tol) ||(iteration == max_accept_iter)

            % Calculate the reactions
            if TangentID == 1
                f_ext_E = K_E * dofs_E + K_EF * dofs_F;
            elseif TangentID == 2
                f_ext_E = DRdu_E * dofs_E + DRdu_EF * dofs_F;
            else
                disp("Check your choice of solver!")
            end
            
            f_external(listofnodes_ebc) = f_ext_E;
            f_external(listofnodes_nbc) = f_ext_F;

            % Reactions
            Reactions_y = sum(abs(f_ext_E(fixnodes_applied(2,:) == 2)))/2;
            Sum_Fx = sum(f_ext_E(fixnodes_applied(2,:) == 1));
            Sum_Fy = sum(f_ext_E(fixnodes_applied(2,:) == 2));

%             disp("End of load increment no.: " + num2str(increment))
%             disp("iteration " + num2str(iteration) + " Res_F_norm: " + Res_F_F_norm_first + "  " + Res_F_F_norm(iteration+1))
%             disp("iteration " + num2str(iteration) + " Res_u_norm: " + Res_u_norm_first   + "  " + Res_u_norm(iteration))
%             disp("========================================")

            break
        else
            iteration = iteration + 1;
        end

    end

elseif IsProj == 1
    % ---------------------------------------------------------------------
    %%%%%%%%%%%%%%%%%%%%% PLOTTING PROPERTIES ROUTINE %%%%%%%%%%%%%%%%%%%%%
    % ---------------------------------------------------------------------

    % Create empty entries for solver variables
    Res_F_F_norm = [];
    Res_u_norm   = [];
    iteration    = [];
    Reactions_y  = [];
    lr_vec_inclist = [];

    for i = 1:ndomains
        % Evaluate stiffness matrix [K] and history variable matrix at the start of the increment
        Domain_Vector(i).history_var_mat_previousinc=Domain_Vector(i).history_var_mat;
        if SolverID == 3 && Domain_Vector(i).DomainDamage
            Domain_Vector(i).damage_mat_inl_prev_inc = Domain_Vector(i).damage_mat_inl_prev_inc_stored;
        end
        Domain_Vector(i) = func_domain_global_stiffness(Domain_Vector(i),IsProj);
    end
else
    disp("Check your IsProj variable")
end

end


