function[OldDomSplit]=func_splitdomain_od(OldDomSplit,Base_Domain,SolverID,j_store, j_not_store, a_vec1, a_vec2)

%create a new domain vector
%assign the elements with no damage to the first position of the new domain
%vector and elements with damage to the second position


OldDomSplit(2).local_strain_mat_stored      = Base_Domain.local_strain_mat_stored(j_store,:);
OldDomSplit(2).nonlocal_strain_mat_stored   = Base_Domain.nonlocal_strain_mat_stored(j_store,:);
OldDomSplit(2).history_var_mat              = Base_Domain.history_var_mat(j_store,:);
OldDomSplit(2).history_var_mat_previousinc  = Base_Domain.history_var_mat_previousinc(j_store,:);

if SolverID == 3
    OldDomSplit(2).damage_mat_inl                   = Base_Domain.damage_mat_inl(j_store,:);
    OldDomSplit(2).damage_mat_inl_prev_inc          = Base_Domain.damage_mat_inl_prev_inc(j_store,:);
    OldDomSplit(2).damage_mat_inl_prev_inc_stored   = Base_Domain.damage_mat_inl_prev_inc_stored(j_store,:);
end



OldDomSplit(2).history_var_mat_stored   = Base_Domain.history_var_mat_stored(j_store,:);
OldDomSplit(2).stress_s1_mat_stored     = Base_Domain.stress_s1_mat_stored(j_store,:);
OldDomSplit(2).damage_mat               = Base_Domain.damage_mat(j_store,:);
OldDomSplit(2).damage_mat_stored        = Base_Domain.damage_mat_stored(j_store,:);
OldDomSplit(2).prop_mat                 = Base_Domain.prop_mat(j_store,:);
OldDomSplit(2).prop_mat_stored          = Base_Domain.prop_mat_stored(j_store,:);


%modify domain 1 to only have properties of the elements remaining
OldDomSplit(1).local_strain_mat_stored      = Base_Domain.local_strain_mat_stored(j_not_store,:);
OldDomSplit(1).nonlocal_strain_mat_stored   = Base_Domain.nonlocal_strain_mat_stored(j_not_store,:);
OldDomSplit(1).history_var_mat              = Base_Domain.history_var_mat(j_not_store,:);
OldDomSplit(1).history_var_mat_previousinc  = Base_Domain.history_var_mat_previousinc(j_not_store,:);


if SolverID == 3
    OldDomSplit(1).damage_mat_inl                   = Base_Domain.damage_mat_inl(j_not_store,:);
    OldDomSplit(1).damage_mat_inl_prev_inc          = Base_Domain.damage_mat_inl_prev_inc(j_not_store,:);
    OldDomSplit(1).damage_mat_inl_prev_inc_stored   = Base_Domain.damage_mat_inl_prev_inc_stored(j_not_store,:);
end


OldDomSplit(1).history_var_mat_stored   = Base_Domain.history_var_mat_stored(j_not_store,:);
OldDomSplit(1).stress_s1_mat_stored     = Base_Domain.stress_s1_mat_stored(j_not_store,:);
OldDomSplit(1).damage_mat               = Base_Domain.damage_mat(j_not_store,:);
OldDomSplit(1).damage_mat_stored        = Base_Domain.damage_mat_stored(j_not_store,:);
OldDomSplit(1).prop_mat                 = Base_Domain.prop_mat(j_not_store,:);
OldDomSplit(1).prop_mat_stored          = Base_Domain.prop_mat_stored(j_not_store,:);



OldDomSplit(2).Res_F=Base_Domain.Res_F(a_vec2);
OldDomSplit(2).f_external=Base_Domain.f_external(a_vec2);
OldDomSplit(2).K=Base_Domain.K(a_vec2,a_vec2);
OldDomSplit(2).DRdu = []; % Base_Domain.DRdu(a_vec2,a_vec2);
OldDomSplit(2).dofs=Base_Domain.dofs(a_vec2);
OldDomSplit(2).dofs_stored=Base_Domain.dofs_stored(a_vec2);


OldDomSplit(1).Res_F=Base_Domain.Res_F(a_vec1);
OldDomSplit(1).f_external=Base_Domain.f_external(a_vec1);
OldDomSplit(1).K=Base_Domain.K(a_vec1,a_vec1);
OldDomSplit(1).DRdu = []; % =Base_Domain.DRdu(a_vec1,a_vec1);
OldDomSplit(1).dofs=Base_Domain.dofs(a_vec1);
OldDomSplit(1).dofs_stored=Base_Domain.dofs_stored(a_vec1);


end