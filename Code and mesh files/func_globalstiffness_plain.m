function [K,DRdu,Res_F,history_var_mat,damage_mat_it,gausspoints_prop_mat,nodes_prop_mat] = func_globalstiffness_plain(domain,DomainDamage,dofs,ndof,nnodes,maxnodes,nelnodes,coords,ncoord,nelem,connect,Delastic,history_var_mat_previousinc,alpha_val,beta_val,e_delta,dmax,n_hood,weights,numelem4node,IsProj)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====================== ASSEMBLE THE TANGENT MATRIX ======================
% ===================== ASSEMBLE THE RESIDUAL VECTOR ======================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include global variables
func_include_flags;

% Calculate the local damage matrix for the nonlocal integral method
if SolverID == 3 && DomainDamage == 1
    damage_mat_it = func_localdamagemat(domain); 
else
    damage_mat_it = [];
end

% tic
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOLVER ROUTINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
if IsProj == 0

    % Create empty entries for element/nodal properties 
    gausspoints_prop_mat = [];
    nodes_prop_mat = [];   
   
    % Initializing global stiffness, global internal force and gather matrices 
    % at zero:
    % K                   = zeros(ndof*nnodes,ndof*nnodes);
    K                   = sparse(ndof*nnodes,ndof*nnodes);

    DRdu                = []; 
    Res_F               = zeros(ndof*nnodes,1);
    residual_elpos      = zeros(maxnodes*ndof,maxnodes*ndof);
    residual_elneg      = zeros(maxnodes*ndof,maxnodes*ndof);
    lmncoord            = zeros(ncoord,maxnodes);
    lmndof              = zeros(ndof,maxnodes);
    history_var_mat     = zeros(nelem,4);




    % Loop over all the elements
    for lmn = 1:nelem
    
        % Extract coords of nodes, DOF for the current element
        for a = 1:nelnodes(lmn)
            for i = 1:ncoord
                lmncoord(i,a) = coords(i,connect(a,lmn));
            end
            for i = 1:ndof
                lmndof(i,a) = dofs(ndof*(connect(a,lmn)-1)+i);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANALYTICAL CALCULATION - TANGENT MATRIX (for each element): 
        if DomainDamage == 0
            [k_el, Res_F_el, history_var_mat(lmn,:), ~, ~] = func_elstif_NoDamage(lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,1,IsProj);
        else
                
            if TangentID == 1
    
                % -------------------------------------------------------------
                if SolverID == 1 
                    [k_el, Res_F_el, history_var_mat(lmn,:), ~, ~] = func_elstif_Local(lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,1,IsProj);
                % -------------------------------------------------------------
                elseif SolverID == 2
                    [k_el, Res_F_el, history_var_mat(lmn,:), ~, ~] = func_elstif_Nonlocgradient(lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),g,alpha_val,beta_val,e_delta,dmax,n_hood,weights,1,IsProj);
                % -------------------------------------------------------------
                elseif SolverID == 3
                    [k_el, Res_F_el, history_var_mat(lmn,:), ~, ~] = func_elstif_Nonlocintegral(lmn,lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,damage_mat_it,1,IsProj);
                else
                    disp("Check your SolverID - globalstiffness")
                end
    
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % NUMERICAL APPROXIMATION - TANGENT MATRIX (for each element):
        
            if TangentID == 2
                    
                % -------------------------------------------------------------
                if SolverID == 1
                    [~, Res_F_el, history_var_mat(lmn,:), ~, ~] = func_elstif_Local(lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,1,IsProj);    
                % -------------------------------------------------------------
                elseif SolverID == 2
                    [~, Res_F_el, history_var_mat(lmn,:), ~, ~] = func_elstif_Nonlocgradient(lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),g,alpha_val,beta_val,e_delta,dmax,n_hood,weights,1,IsProj);
                % -------------------------------------------------------------
                elseif SolverID == 3
                    [~, Res_F_el, history_var_mat(lmn,:), ~, ~] = func_elstif_Nonlocintegral(lmn,lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,damage_mat_it,1,IsProj);
                % -------------------------------------------------------------
                else 
                    disp("Check your SolverID - globalstiffness")
                end
    
                % element stiffness (based on residual at plus/minus lmndof - central difference)
                dl = 1e-10; % Tolerance 
                col_id = 0; % Counter of columns in dRdu (range: 1-12)
            
                for a = 1:nelnodes(lmn)
                    for i = 1:ndof
                        
                        % Increase counter
                        col_id = col_id + 1;
            
                        % Compute residual_el at (u+dl(dof)) of size: 12x1, and store 
                        % it in the next column of residual
                        lmndofi = lmndof;    
                        lmndofi(i,a) = lmndofi(i,a) + dl;
                        % -----------------------------------------------------
                        if SolverID == 1
                            [~, residual_elpos(:,col_id), ~, ~, ~] = func_elstif_Local(lmncoord,lmndofi,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,0,IsProj);
                        % -----------------------------------------------------
                        elseif SolverID == 2
                            [~, residual_elpos(:,col_id), ~, ~, ~] = func_elstif_Nonlocgradient(lmncoord,lmndofi,Delastic,history_var_mat_previousinc(lmn,:),g,alpha_val,beta_val,e_delta,dmax,n_hood,weights,0,IsProj); % Final size: 12x12
                        % -----------------------------------------------------
                        elseif SolverID == 3
                            [~, residual_elpos(:,col_id), ~, ~, ~] = func_elstif_Nonlocintegral(lmn,lmncoord,lmndofi,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,damage_mat_it,0,IsProj);
                        % -----------------------------------------------------
                        else
                            disp("Check your SolverID - globalstiffness")
                        end
    
                        % Compute residual_el at (u-dl(dof)) of size: 12x1, and store 
                        % it in the next column of residual
                        lmndofi = lmndof;    
                        lmndofi(i,a) = lmndofi(i,a) - dl;
                        % -----------------------------------------------------
                        if SolverID == 1 
                            [~, residual_elneg(:,col_id), ~, ~, ~] = func_elstif_Local(lmncoord,lmndofi,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,0,IsProj);
                        % -----------------------------------------------------
                        elseif SolverID == 2
                            [~, residual_elneg(:,col_id), ~, ~, ~] = func_elstif_Nonlocgradient(lmncoord,lmndofi,Delastic,history_var_mat_previousinc(lmn,:),g,alpha_val,beta_val,e_delta,dmax,n_hood,weights,0,IsProj); % Final size: 12x12
                        % -----------------------------------------------------
                        elseif SolverID == 3
                            [~, residual_elneg(:,col_id), ~, ~, ~] = func_elstif_Nonlocintegral(lmn,lmncoord,lmndofi,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,damage_mat_it,0,IsProj);
                        % -----------------------------------------------------
                        else
                            disp("Check your SolverID - globalstiffness")
                        end
                    end
                end
     
                % Compute "partial R over partial u" for all dofs at one step 
                dRdu_el = 1/(2*dl) * (residual_elpos - residual_elneg); % Size: 12x12
    
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ASSEMBLE GLOBAL MATRICES: 
        % a) analytical tangent K matrix
        % b) numerical tangent DRdu matrix
        % c) residual vector
        
        
        % a) Analytical tangent matrix
        if DomainDamage==0
            for a = 1:nelnodes(lmn)
                for i = 1:ndof
                    for b = 1:nelnodes(lmn)
                        for k = 1:ndof
                            rw = ndof*(connect(a,lmn)-1)+i;
                            cl = ndof*(connect(b,lmn)-1)+k;
                            K(rw,cl) = K(rw,cl) + k_el(ndof*(a-1)+i,ndof*(b-1)+k);
                        end
                    end
                end
            end
        else
            if TangentID == 1
                for a = 1:nelnodes(lmn)
                    for i = 1:ndof
                        for b = 1:nelnodes(lmn)
                            for k = 1:ndof
                                rw = ndof*(connect(a,lmn)-1)+i;
                                cl = ndof*(connect(b,lmn)-1)+k;
                                K(rw,cl) = K(rw,cl) + k_el(ndof*(a-1)+i,ndof*(b-1)+k);
                            end
                        end
                    end
                end
            end
        
%             % b) Numerical tangent matrix
%             if TangentID  == 2
%                 for a = 1:nelnodes(lmn)
%                     for i = 1:ndof
%                         for b = 1:nelnodes(lmn)
%                             for k = 1:ndof
%                                 rw = ndof*(connect(a,lmn)-1)+i;
%                                 cl = ndof*(connect(b,lmn)-1)+k;
%                                 DRdu(rw,cl) = DRdu(rw,cl) + dRdu_el(ndof*(a-1)+i,ndof*(b-1)+k);
%                             end
%                         end
%                     end
%                 end
%             end
        end
        % c) Residual vector
        for a = 1:nelnodes(lmn)
            for i = 1:ndof
                rw = ndof*(connect(a,lmn)-1)+i;
                Res_F(rw,1) = Res_F(rw,1) + Res_F_el(ndof*(a-1)+i,1);
            end
        end
    
    end

% toc

elseif IsProj == 1
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%% PLOTTING PROPERTIES ROUTINE %%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------

    % Create empty entries for solver variables
    K         = [];
    DRdu      = [];
    Res_F     = [];
    history_var_mat = [];

    % Initializing element/nodal properties and gather matrices at zero
    gausspoints_prop_mat = zeros(nelem * func_numberofintegrationpoints,9);
    nodes_prop_mat       = zeros(nnodes,9);
    lmncoord             = zeros(ncoord,maxnodes);
    lmndof               = zeros(ndof,maxnodes);


    % Loop over all the elements
    for lmn = 1:nelem
    
        % Extract coords of nodes, DOF for the current element
        for a = 1:nelnodes(lmn)
            for i = 1:ncoord
                lmncoord(i,a) = coords(i,connect(a,lmn));
            end
            for i = 1:ndof
                lmndof(i,a) = dofs(ndof*(connect(a,lmn)-1)+i);
            end
        end

        if DomainDamage == 0
            [~,~,~,gausspoints_prop_mat(4*lmn-3:4*lmn,:), nodes_prop_mat_elem] = func_elstif_NoDamage(lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,1,IsProj);
       
        else


            % Compute element/nodal properties
            % -----------------------------------------------------------------
            if SolverID == 1
                [~, ~, ~, gausspoints_prop_mat(4*lmn-3:4*lmn,:), nodes_prop_mat_elem] = func_elstif_Local(lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,1,IsProj);
            % -----------------------------------------------------------------
            elseif SolverID == 2
                [~, ~, ~, gausspoints_prop_mat(4*lmn-3:4*lmn,:), nodes_prop_mat_elem] = func_elstif_Nonlocgradient(lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),g,alpha_val,beta_val,e_delta,dmax,n_hood,weights,1,IsProj);
            % -----------------------------------------------------------------
            elseif SolverID == 3
                [~, ~, ~, gausspoints_prop_mat(4*lmn-3:4*lmn,:), nodes_prop_mat_elem] = func_elstif_Nonlocintegral(lmn,lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,damage_mat_it,1,IsProj);
            % -----------------------------------------------------------------
            else
                disp("Check your SolverID - globalstiffness")
            end
        end
        
        % Populate the properties matrix
        for a = 1:nelnodes(lmn)
            rw = connect(a,lmn);
            nodes_prop_mat(rw,:) = nodes_prop_mat(rw,:) + nodes_prop_mat_elem(a,:);
        end

    end

    % Calculate the final matrix of nodal properties
    nodes_prop_mat = nodes_prop_mat ./ numelem4node;

else 

    disp("Check your IsProj variable - globalstiffness")

end


end
