%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start another kernel of matlab and share its engine 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
format compact
format long

% Sample file used to connect the engine (provided)
name_of_file = "1dom2.txt";
modeOfFailure = "1 2 0 5";

% pyrunfile(pwd + "/python_files/read_input_file.py", filename = name_of_file, path = pwd, modeOfFailure = modeOfFailure);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======================= INCLUDE GLOBAL VARIABLES ========================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
func_include_flags;
load custom_colormaps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ========================== INITIAL PARAMETERS ===========================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MATLAB PARAMETERS
newDecomp       = 1; % This will be used by the schur complement to know whether it must recalculate inverses in some domain
adadode         = 0; % Do adaptive domain decomp (1) or not (0)
j_store         = []; % Identities of elements in damage zone
j_not_store     = []; % Identities of elements not in damage zone
a_vec1          = []; % Identities of dofs not in damage zone
a_vec2          = []; % Identities of dofs in damage zone
converged       = 0; % Variable for checking whether a time step has converged
ALC             = 0; % Use arclength control or not
plot_boundaries = 1; % Plot location of decomposition boundary

% PYTHON PARAMETERS
first_time      = 1; % Used by python to determine whether decomposition has started
elastic_inc     = 1; % Used by python to identify the redundant border in the elastic zone
wx              = [1,1]; % Weight of distance in the x direction. used by python
wy              = [1,1]; % Weight of distance in the y direction. used by python
num_cont        = 0; % Number of contours detected
bound_damage    = [0,0]; % Variable initializer for what previous damage contours
centroids       = [0,0]; % Centroids of the contours detected
orrcont         = []; % Real coordinates of contours created
Khh_inv         = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====================== CREATE INPUT FILE AND LOAD IT=====================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_name = "ADNT_Coarse"; % Optional model name

infile = fopen(strcat(model_name, "_parameters.txt"),'r');
[SolverID,TangentID,ncoord,ndof,lc, ...
    increment,inc_success_counter,min_iter,max_iter,max_accept_iter, ...
    loadfactor,dlfactor,dlfactor_incr_threshold,increment_plot_threshold,loadfactor_plot_threshold, ...
    flaglf,countflaglf,incrflag,flagplot, ...
    ndomains,Domain_Vector] = func_read_input_file(infile,model_name);
fclose(infile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loadfactor_plot_threshold   = 0; %This was set to zero to allow for comparison between 1 domain and multiple domains
Base_Domain                 = Domain_Vector; % Base domain which serves as reference before domains are split
Coords                      = Base_Domain.coords;
[K_pos_mat, fixnodes]       = func_bound_pos(Domain_Vector,ndomains); % Positions of boundary of various domains
undamaged_elements_stored   = Base_Domain.nelem;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ================== SPECIFY VALUES FOR GLOBAL VARIABLES ==================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Non-local gradient parameter (g = lc^2/2)
g = lc^2/2;

% -------------------------------------------------------------------------
% Choose solution scheme
% SolverID: 1 - Local
%           2 - Nonlocal Gradient
%           3 - Nonlocal Integral
% TangentID: 1 - Analytical
%            2 - Numerical

% -------------------------------------------------------------------------
if SolverID == 1
    solver_string = "Local";
elseif SolverID == 2
    solver_string = "Nonlocal_gradient";
elseif SolverID == 3
    solver_string = "Nonlocal_integral";
else
    disp("Check your SolverID - FEM Main Script")
end

if TangentID == 1
    tangent_string = "Analytical";
elseif TangentID == 2
    tangent_string = "Numerical";
else
    disp("Check your TangentID - FEM_Main_Script")
end

% -------------------------------------------------------------------------
% Calculate matrix of neighbouring points and their weights for the nonlocal integral method
if SolverID == 3
    for i = 1:ndomains
        if (Domain_Vector(i).DomainDamage == 1)
            [Domain_Vector(i).n_hood, Domain_Vector(i).weights] = func_nhood_gausspts(lc, Domain_Vector(i));
        end
    end
end

lr      = 0; % Damping parameter
flag_lr = 0; % Ensures damping parameter isn't stuck at its initial value
converge_counter = 0; % Used to check whether to increase the damping parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ===================== NEWTON - RAPHSON ANALYSIS =====================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NR_time = 0; print_time = 0; dd_time = 0; load_flag = 0; counter_any_increment = 0;

while loadfactor <= 1.01 && countflaglf == 0 && incrflag < 100

    % Condition to terminate the while loop: the total load is applied
    if flaglf == true; countflaglf = 1; end

    % Condition to terminate the while loop: dlfactor is too small
    if dlfactor < 10^-40; break; end

    % Comment out the next block of code to do a purely Newton-Raphson analysis.
    if dlfactor < 10^-4  && flag_lr == 0
        lr = 1E-9;
        flag_lr = 1;
    end
    
    % ---------------------------------------------------------------------
    % --------------------------- SOLVER ROUTINE --------------------------
    % ---------------------------------------------------------------------

    IsProj = 0;

    % Perform the Newton Raphson analysis
    disp("========================================")
    disp(strcat("Loadfactor: " + loadfactor))
    disp(strcat("dlfactor: " + dlfactor))
    disp(strcat("lr: " + lr))    

    % This is used purely for plotting purposes
    Domain_Vector_before = Domain_Vector; 
    
    % ---------------------------------------------------------------------
    % This counter increases at each increment, no matter what (whether NR
    % fails to converge or if dd is repeated)
    counter_any_increment = counter_any_increment + 1;
    % ---------------------------------------------------------------------

    start_NR = tic; % Time at start of NR routine
    [Domain_Vector, Res_F_F_norm, Res_u_norm, last_iteration, Reactions_y,Khh_inv,lr_vec_inclist,Beta] = func_NewtonRaphson_DispControl([fixnodes(1:2,:); fixnodes(3,:)*loadfactor],increment,Domain_Vector,K_pos_mat,ndomains,IsProj,lr,Khh_inv);
    
    if length(Domain_Vector) == 2
        max(Domain_Vector_graph(2).gausspoints_prop_mat(:,1))    
    end

    % Total time taken by NR routine
    NR_time_any_increment(counter_any_increment)        = toc(start_NR);
    iterations_at_any_increment(counter_any_increment)  = last_iteration;
    dlfactor_vec_any_increment(counter_any_increment)   = dlfactor;

    NR_time_increment(increment) = toc(start_NR);
    NR_time                      = NR_time + toc(start_NR);

    disp("iteration " + num2str(last_iteration))

    % Store Reactions
    if last_iteration < max_accept_iter
        Reactions_mat(inc_success_counter,1)    = loadfactor;
        Reactions_mat(inc_success_counter,2)    = Reactions_y;
        lr_vec(inc_success_counter)             = lr;
        dlfactor_vec(inc_success_counter)       = dlfactor;
        Res_u{inc_success_counter}              = Res_u_norm;
        Res_F_F{inc_success_counter}            = Res_F_F_norm;
        lr_vec_inclist_success_counter{inc_success_counter} = lr_vec_inclist;
        Beta_vec{inc_success_counter}           = Beta;
    end

    if inc_success_counter > 1
        Domain_Vector_icr       = Domain_Vector;
        flaglf_icr              = flaglf;
        incrflag_icr            = incrflag;
        loadfactor_icr          = loadfactor_stored;
        dlfactor_icr            = dlfactor;
        inc_success_counter_icr = inc_success_counter;
    end

    % -----------------------------------------------------------------
    % ----------------- PLOTTING PROPERTIES ROUTINE -------------------
    % -----------------------------------------------------------------

    % Time at start of printing routine
    disp("Printing routine")
    start_print_contours = tic; 

    if (last_iteration ~= max_accept_iter) && (loadfactor > loadfactor_plot_threshold) % (rem(increment,increment_plot_threshold) == 0)

        IsProj = 1;

        % Save plots to be used for domain decomp
        if mod(increment,1) == 0 || adadode == 1
            [Domain_Vector_graph, ~,~,~,~,~,~,~] = func_NewtonRaphson_DispControl([fixnodes(1:2,:); fixnodes(3,:)*loadfactor],increment,Domain_Vector,K_pos_mat,ndomains,IsProj,lr,Khh_inv);
            func_saveplots(ndomains,model_name,inc_success_counter,last_iteration,increment,loadfactor,Res_F_F_norm,Res_u_norm,solver_string,tangent_string,Domain_Vector_graph);
%             max_damage = max(Domain_Vector_graph.gausspoints_prop_mat(:,1))
        end
    end

    % Total time taken by printing routine
    print_time_increment(increment) = toc(start_print_contours); 
    print_time = print_time + toc(start_print_contours); 


    % ---------------------------------------------------------------------
    % ------------------------ ADAPT LOAD ROUTINE -------------------------
    % ---------------------------------------------------------------------
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if loadfactor >= 0.4 && load_flag == 0
        dlfactor = 0.01;
        dlfactor_incr_threshold = 0.01;
        load_flag = 1;      
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Adapt load incrementation w.r.t. the number of iterations needed for convergence
    if last_iteration <= min_iter
        % Fast Convergence - dlfactor increases
        [Domain_Vector,loadfactor_stored,incrflag,flagplot,loadfactor,increment,inc_success_counter,flaglf,dlfactor]...
            = func_adapt_quick_convergence(ndomains,Domain_Vector,loadfactor,last_iteration,dlfactor,increment,inc_success_counter);
        converged = 1;
        converge_counter = 0;

    elseif (min_iter < last_iteration) && (last_iteration <= max_iter)
        % Moderate Convergence - dlfactor remains the same
        [Domain_Vector,loadfactor_stored,incrflag,flagplot,loadfactor,increment,inc_success_counter,flaglf]...
            = func_adapt_moderate_convergence(ndomains,Domain_Vector,loadfactor,dlfactor,increment,inc_success_counter);
        converged = 1;
        converge_counter = 0;

    elseif (max_iter < last_iteration) && (last_iteration < max_accept_iter)
        % Slow Convergence - dlfactor decreases
        [Domain_Vector,loadfactor_stored,incrflag,flagplot,loadfactor,increment,inc_success_counter,flaglf,dlfactor] ...
            = func_adapt_slow_convergence(ndomains,Domain_Vector,loadfactor,last_iteration,dlfactor,increment,inc_success_counter);
        converged = 1;
        converge_counter = 0;

    else
        % No Convergence - discard the last step and repeat with smaller load value
        disp("what happened here")
        [Domain_Vector,flagplot,flaglf,countflaglf,incrflag,loadfactor,dlfactor]...
            = func_adapt_no_convergence(ndomains,Domain_Vector,incrflag,loadfactor_stored,dlfactor);
        % disp(strcat("loadfactor" + loadfactor))
        % disp(strcat("dlfactor" + dlfactor))
        % disp(strcat("lr" + lr))
        converged = 0;
        converge_counter = converge_counter + 1;

        if lr > 0 && converge_counter > 5
            % if the problem has not converged for 5 consecutive times, increase the damping parameter
            lr = min(lr * 10, 1e-3); % lr * 10
            converge_counter = 0;
        end
    end
    
    if elastic_inc == 1
        coordinates_redundant_border_stored = [];
    else
        coordinates_redundant_border_stored = coordinates_redundant_border;
    end

    if adadode == 1 && converged == 1
        disp("Image segmentation and domain decomposition routine")
        start_dd = tic; % Time at start of DD routine

        % Check for contours and extract information
        pure_python_dd = tic;

        [repeat,pts_in_damage,bound_damage,decompornot,wx,wy,num_cont,first_time,elastic_inc,orrcont,centroids,coordinates_redundant_border]...
            = pyrunfile(pwd + "/python_files/find_contours.py",["repeat","output","new_decomp",...
            "new_decompbool","weightx","weighty","num_of_conts","first_time","elastic_inc",...
            "contours_send_back","new_centroid","coordinates_redundant_border","increment_py"], path = strcat(pwd+"/" ...
            + model_name + "_DamContour_" + solver_string + ...
            "_inc_" + int2str(increment-1) + ".png"), ... 
            coords                       = py.numpy.array(Coords),...
            prev_number_of_contours      = num_cont, ... 
            old_small_cnt                = py.numpy.array(bound_damage),...
            LC                           = lc, ...
            weightx                      = py.numpy.array(wx), ... 
            weighty                      = py.numpy.array(wy),...
            first_time                   = first_time, ... 
            elastic_inc                  = elastic_inc, ...
            centroids                    = py.numpy.array(centroids), ...
            coordinates_redundant_border_stored = coordinates_redundant_border_stored,...
            increment_py                 = increment - 1);

            dd_time_pure_python_any_increment(counter_any_increment) = toc(pure_python_dd);

            if decompornot % decomp to 2 domains or not
                pts_in_damage = cellfun(@double,cell(pts_in_damage));

                if repeat == 1
                    Domain_Vector = Domain_Vector_icr;
                end
    
                if ~isempty(pts_in_damage) % Only do the decomposition if there are points inside the damaged region
%                     disp("Here 1")
                    [Domain_Vector,Base_Domain,ndomains,K_pos_mat,fixnodes,j_store, j_not_store,a_vec1, a_vec2] = func_moveElements(ndomains,pts_in_damage,Domain_Vector,Base_Domain,j_store, j_not_store, a_vec1, a_vec2,lc,SolverID);
                end

            else
%                 disp("Here 2")
                if first_time == 0
%                     disp("Here 3")
                    [Base_Domain,~] = func_combinedomains(ndomains,Base_Domain,Domain_Vector,j_store, j_not_store, a_vec1, a_vec2);
                    [Domain_Vector] = func_splitdomain_od(Domain_Vector,Base_Domain,SolverID,j_store, j_not_store, a_vec1, a_vec2);
                end
            end
            
            if undamaged_elements_stored ~= Domain_Vector(1).nelem
                newDecomp = 1;
                undamaged_elements_stored = Domain_Vector(1).nelem;
            else
                newDecomp = 0;
            end
            
            if length(Domain_Vector) == 2
                damaged_nnodes_any_increment(counter_any_increment) = Domain_Vector(2).nnodes;
            end

        % Total time taken by DD routine
        dd_time_increment(increment) = toc(start_dd);
        dd_time_any_increment(counter_any_increment) = toc(start_dd);
        dd_time = dd_time + toc(start_dd);
        
        % If damage crossed the old damage interface, redo the loadstep
        if repeat == 1
            increment           = increment - 1;
            incrflag            = incrflag_icr;
            loadfactor_stored   = loadfactor_icr;
            dlfactor            = dlfactor_icr;
            inc_success_counter = inc_success_counter_icr;

            for i = 1:ndomains
                Domain_Vector(i).dofs = Domain_Vector(i).dofs_stored;
                Domain_Vector(i).history_var_mat = Domain_Vector(i).history_var_mat_stored;
                if SolverID == 3 && Domain_Vector(i).DomainDamage
                    Domain_Vector(i).damage_mat_inl = Domain_Vector(i).damage_mat_inl_prev_inc_stored;
                end
            end

            % Allow the while loop to continue running in case loadfactor has hit the value 1
            flagplot    = 0;
            flaglf      = false;
            countflaglf = 0;
            incrflag    = incrflag + 1;

            % Use stored values for load factor (kappa and damage are already being stored)
            loadfactor = loadfactor_stored;

            % Reduce the applied load
            loadfactor = min(loadfactor + dlfactor,1);

            % No Convergence - discard the last step and repeat with smaller load value
            converged = 0;

        end
        
        % Show the decomposition boundary if plot_boundaries is 1
        if first_time == 0 && plot_boundaries == 1 && converged == 1

            [Domain_Vector_graph, ~,~,~,~,~,~] = func_NewtonRaphson_DispControl([fixnodes(1:2,:); fixnodes(3,:)*loadfactor],increment-1,Domain_Vector,K_pos_mat,ndomains,1,lr,Khh_inv);
            figure('visible','off'); axis equal; axis off; caxis([0,1]); hold on; 
            colormap(jet);  % https://www.mathworks.com/help/matlab/ref/colormap.html

            for i = 1:ndomains
                Domain_Vector_graph(i).plotcontours(model_name,inc_success_counter,last_iteration,increment,loadfactor,Res_F_F_norm,Res_u_norm,solver_string,tangent_string,'k')
            end

            % Plot scatter of unhealthy nodes
            % scatter(Domain_Vector(2).coords(1,:),Domain_Vector(2).coords(2,:),14,'x','MarkerEdgeColor',[0.9 0.7 0.1]);
            % hold on

            % Rearranging nodes in current element to form closed quadrilateral
            forcoordx = Domain_Vector(1).coords(1,K_pos_mat(1,K_pos_mat(1,:) ~= 0));
            forcoordy = Domain_Vector(1).coords(2,K_pos_mat(1,K_pos_mat(1,:) ~= 0));
            scatter(forcoordx,forcoordy,'.w')

            origicont = cellfun(@double,cell(orrcont),UniformOutput = false);       
            for i = 1:size(origicont,2)
                plotorr = origicont(i);
                plotorr = plotorr{1};
                plot([plotorr(:,:,1);plotorr(1,:,1)],[plotorr(:,:,2);plotorr(1,:,2)],'r')
                hold on
            end

%             % Plot bound_damage
%             double_damage_double = double(bound_damage);
%             plot([double_damage_double(1,:,1) double_damage_double(1,1,1)],[double_damage_double(1,:,2) double_damage_double(1,1,2)],'-.m','LineWidth',1.6);
%             plot([double_damage_double(2,:,1) double_damage_double(2,1,1)],[double_damage_double(2,:,2) double_damage_double(2,2,2)],'-.m','LineWidth',1.6);

%             xlim([20 42]); ylim([15 36]);

            name2 = strcat(pwd + "\" + model_name + "_boundary_" + solver_string + "_inc_" + int2str(increment-1) + ".png");
            exportgraphics(gcf,name2,'Resolution',300)
            hold off
            close all;
        end
    end
end


fprintf("Time(s) spent for just NR: %.4f\n\n", NR_time)
fprintf("Time(s) spent for just dd: %.4f\n\n", dd_time)
fprintf("Time(s) spent for printing images: %.4f\n\n", print_time)
fprintf("Time(s) spent for both NR and dd combined: %.4f\n\n", NR_time + dd_time)
fprintf("Time(s) spent for NR, dd and printing images combined: %.4f\n\n", NR_time + dd_time + print_time)

% -------------------------------------------------------------------------
% Plotting reactions
Reactions_mat = [0 0; Reactions_mat];

if adadode == 0
    dd_time_increment = 0;
end
save Results_ADNT_A1_final_tol1e5.mat Domain_Vector_graph NR_time dd_time print_time Reactions_mat Res_u Res_F_F NR_time_increment print_time_increment dd_time_increment lr_vec dlfactor_vec NR_time_any_increment iterations_at_any_increment lr_vec_inclist_success_counter Beta_vec dd_time_any_increment damaged_nnodes_any_increment dd_time_pure_python_any_increment dlfactor_vec_any_increment



