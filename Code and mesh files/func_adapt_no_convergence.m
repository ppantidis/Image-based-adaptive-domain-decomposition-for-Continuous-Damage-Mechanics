function [Domain_Vector,flagplot,flaglf,countflaglf,incrflag,loadfactor,dlfactor] = func_adapt_no_convergence(ndomains,Domain_Vector,incrflag,loadfactor_stored,dlfactor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======== ADAPT LOADFACTOR WHEN NR CONVERGENCE IS NOT CONVERGING =========
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include global variables
func_include_flags;

for i=1:ndomains
            Domain_Vector(i).dofs=Domain_Vector(i).dofs_stored;
            Domain_Vector(i).history_var_mat=Domain_Vector(i).history_var_mat_stored;
            if SolverID == 3 && Domain_Vector(i).DomainDamage
                Domain_Vector(i).damage_mat_inl=Domain_Vector(i).damage_mat_inl_prev_inc_stored;
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
dlfactor = dlfactor/2;
loadfactor = min(loadfactor + dlfactor,1);

disp("you just entered highly nonlinear, check results")
disp("========================================")        

end