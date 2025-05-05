function[Domain_Vector,Base_Domain,ndomains,K_pos,fixnodes,j_store, j_not_store, a_vec1, a_vec2]=func_moveElements(ndomains,pts_in_damage,Domain_Vector,Base_Domain,j_store, j_not_store, a_vec1, a_vec2,lc,SolverID)

% If there is just one domain, then split into two
if ndomains == 1
    Base_Domain = Domain_Vector;
    [Domain_Vector,ndomains,K_pos,fixnodes,j_store,j_not_store,a_vec1,a_vec2] = func_splitdomain(ndomains,pts_in_damage,Base_Domain,lc,SolverID);
else
    % Combine both and then split using new boundary
    [Base_Domain,ndomains]                                                    = func_combinedomains(ndomains,Base_Domain,Domain_Vector, j_store, j_not_store, a_vec1, a_vec2);
    [Domain_Vector,ndomains,K_pos,fixnodes,j_store,j_not_store,a_vec1,a_vec2] = func_splitdomain(ndomains,pts_in_damage,Base_Domain,lc,SolverID);
end

end