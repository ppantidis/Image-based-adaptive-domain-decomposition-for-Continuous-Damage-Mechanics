function [Domain_Vector,loadfactor_stored,incrflag,flagplot,loadfactor,increment, ...
            inc_success_counter,flaglf,dlfactor] = func_adapt_quick_convergence(ndomains,Domain_Vector,loadfactor,last_iteration,dlfactor,increment,inc_success_counter)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== ADAPT MATRICES AND LOADFACTOR WHEN NR CONVERGENCE IS QUICK =======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include global variables
func_include_flags;

% Save variables in case of not convergence of next increment
for i=1:ndomains
    Domain_Vector(i).dofs_stored=Domain_Vector(i).dofs;
    Domain_Vector(i).history_var_mat_stored=Domain_Vector(i).history_var_mat;
    Domain_Vector(i).history_var_mat_stored=Domain_Vector(i).history_var_mat;
    if SolverID == 3 && Domain_Vector(i).DomainDamage
        Domain_Vector(i).damage_mat_inl_prev_inc_stored=Domain_Vector(i).damage_mat_inl;
    end
    
end

loadfactor_stored = loadfactor;


% Re-initialize at 1 the number of attempts for the increment
incrflag = 1; 
flagplot = 1;

% Increment the load
dlfactor            = min(min_iter/last_iteration * dlfactor, dlfactor_incr_threshold);
loadfactor          = min(loadfactor + dlfactor,1);
increment           = increment + 1;
inc_success_counter = inc_success_counter + 1;

% Check if the total load is applied
if loadfactor == 1 
    flaglf = true; 
else
    flaglf = false;
end


end













