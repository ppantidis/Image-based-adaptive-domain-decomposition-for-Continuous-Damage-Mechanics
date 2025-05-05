function[Base_Domain,ndomains]=func_combinedomains(ndomains,Base_Domain,Domain_Vector, j_store, j_not_store, a_vec1, a_vec2)
%this function updates the properties that were split between domains

% Include global variables
func_include_flags;

Base_Domain.damage_mat(j_not_store,:)=Domain_Vector(1).damage_mat;
Base_Domain.damage_mat_stored(j_not_store,:)=Domain_Vector(1).damage_mat_stored;
Base_Domain.local_strain_mat_stored(j_not_store,:)=Domain_Vector(1).local_strain_mat_stored;
Base_Domain.nonlocal_strain_mat_stored(j_not_store,:)=Domain_Vector(1).nonlocal_strain_mat_stored;
Base_Domain.history_var_mat_stored(j_not_store,:)=Domain_Vector(1).history_var_mat_stored;
Base_Domain.history_var_mat_previousinc(j_not_store,:)=Domain_Vector(1).history_var_mat_previousinc;

% if SolverID == 3
%     Base_Domain.damage_mat_inl(j_not_store,:)=Domain_Vector(1).damage_mat_inl;
%     Base_Domain.damage_mat_inl_prev_inc(j_not_store,:)=Domain_Vector(1).damage_mat_inl_prev_inc;
%     Base_Domain.damage_mat_inl_prev_inc_stored(j_not_store,:)=Domain_Vector(1).damage_mat_inl_prev_inc_stored;
% end

Base_Domain.stress_s1_mat_stored(j_not_store,:)=Domain_Vector(1).stress_s1_mat_stored;
Base_Domain.prop_mat(j_not_store,:)=Domain_Vector(1).prop_mat;
Base_Domain.prop_mat_stored(j_not_store,:)=Domain_Vector(1).prop_mat_stored;


Base_Domain.damage_mat(j_store,:)=Domain_Vector(2).damage_mat;
Base_Domain.damage_mat_stored(j_store,:)=Domain_Vector(2).damage_mat_stored;
Base_Domain.local_strain_mat_stored(j_store,:)=Domain_Vector(2).local_strain_mat_stored;
Base_Domain.nonlocal_strain_mat_stored(j_store,:)=Domain_Vector(2).nonlocal_strain_mat_stored;
Base_Domain.history_var_mat_stored(j_store,:)=Domain_Vector(2).history_var_mat_stored;
Base_Domain.history_var_mat_previousinc(j_store,:)=Domain_Vector(2).history_var_mat_previousinc;

if SolverID == 3
    Base_Domain.damage_mat_inl(j_store,:)=Domain_Vector(2).damage_mat_inl;
    Base_Domain.damage_mat_inl_prev_inc(j_store,:)=Domain_Vector(2).damage_mat_inl_prev_inc;
    Base_Domain.damage_mat_inl_prev_inc_stored(j_store,:)=Domain_Vector(2).damage_mat_inl_prev_inc_stored;
end

Base_Domain.stress_s1_mat_stored(j_store,:)=Domain_Vector(2).stress_s1_mat_stored;
Base_Domain.prop_mat(j_store,:)=Domain_Vector(2).prop_mat;
Base_Domain.prop_mat_stored(j_store,:)=Domain_Vector(2).prop_mat_stored;




Base_Domain.Res_F(a_vec1)=Domain_Vector(1).Res_F;
Base_Domain.f_external(a_vec1)=Domain_Vector(1).f_external;
Base_Domain.K(a_vec1,a_vec1)=Domain_Vector(1).K;
Base_Domain.DRdu = []; % (a_vec1,a_vec1)=Domain_Vector(1).DRdu;
Base_Domain.dofs(a_vec1)=Domain_Vector(1).dofs;
Base_Domain.dofs_stored(a_vec1)=Domain_Vector(1).dofs_stored;

Base_Domain.Res_F(a_vec2)=Domain_Vector(2).Res_F;
Base_Domain.f_external(a_vec2)=Domain_Vector(2).f_external;
Base_Domain.K(a_vec2,a_vec2)=Domain_Vector(2).K;
Base_Domain.DRdu = []; % (a_vec2,a_vec2)=Domain_Vector(2).DRdu;
Base_Domain.dofs(a_vec2)=Domain_Vector(2).dofs;
Base_Domain.dofs_stored(a_vec2)=Domain_Vector(2).dofs_stored;

ndomains=ndomains-1;
end