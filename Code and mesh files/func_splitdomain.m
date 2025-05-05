function[NewDomSplit,ndomains,K_pos,fixnodes, j_store, j_not_store, a_vec1, a_vec2]=func_splitdomain(ndomains,pts_in_damage,Base_Domain,lc,SolverID)
%create a new domain vector
%assign the elements with no damage to the first position of the new domain
%vector and elements with damage to the second position

NewDomSplit=domain.empty(ndomains,0);
NewDomSplit(1)=domain;
NewDomSplit(2)=domain;

NewDomSplit(1).DomainDamage=0;
NewDomSplit(2).DomainDamage=1;

%copy common properties
NewDomSplit(1).nprops=Base_Domain.nprops;
NewDomSplit(1).materialprops=Base_Domain.materialprops;
NewDomSplit(1).ncoord=Base_Domain.ncoord;
NewDomSplit(1).ndof=Base_Domain.ndof;
NewDomSplit(1).maxnodes=Base_Domain.maxnodes;
NewDomSplit(1).Delastic=Base_Domain.Delastic;
NewDomSplit(1).alpha_val=Base_Domain.alpha_val;
NewDomSplit(1).beta_val=Base_Domain.beta_val;
NewDomSplit(1).e_delta=Base_Domain.e_delta;
NewDomSplit(1).dmax=Base_Domain.dmax;

NewDomSplit(2).nprops=Base_Domain.nprops;
NewDomSplit(2).materialprops=Base_Domain.materialprops;
NewDomSplit(2).ncoord=Base_Domain.ncoord;
NewDomSplit(2).ndof=Base_Domain.ndof;
NewDomSplit(2).maxnodes=Base_Domain.maxnodes;
NewDomSplit(2).Delastic=Base_Domain.Delastic;
NewDomSplit(2).alpha_val=Base_Domain.alpha_val;
NewDomSplit(2).beta_val=Base_Domain.beta_val;
NewDomSplit(2).e_delta=Base_Domain.e_delta;
NewDomSplit(2).dmax=Base_Domain.dmax;

%create a variable for storing elements to be transferred to domain 2
j_store=[];
size_connect=size(Base_Domain.connect,2);
for i =1:length(pts_in_damage)
    for j =1:size_connect
        if ismember(pts_in_damage(i),Base_Domain.connect(:,j))
            j_store=[j_store,j];
        end
    end
end

%there will be duplicates because of how the code is written. So
%eliminate them
[j_store,~,~]=unique(j_store);

%based on the elements to be transferred, note the elements that are
%going to remain in domain 1
j_not_store=setdiff(Base_Domain.elident_vec, j_store);

%transfer the necessary elements and their properties to domain 2
NewDomSplit(2).connect=Base_Domain.connect(:,j_store);
NewDomSplit(2).nelnodes=Base_Domain.nelnodes(j_store);
NewDomSplit(2).local_strain_mat_stored=Base_Domain.local_strain_mat_stored(j_store,:);
NewDomSplit(2).nonlocal_strain_mat_stored=Base_Domain.nonlocal_strain_mat_stored(j_store,:);
NewDomSplit(2).history_var_mat=Base_Domain.history_var_mat(j_store,:);
NewDomSplit(2).history_var_mat_previousinc=Base_Domain.history_var_mat_previousinc(j_store,:);

if SolverID == 3
    NewDomSplit(2).damage_mat_inl=Base_Domain.damage_mat_inl(j_store,:);
    NewDomSplit(2).damage_mat_inl_prev_inc=Base_Domain.damage_mat_inl_prev_inc(j_store,:);
    NewDomSplit(2).damage_mat_inl_prev_inc_stored=Base_Domain.damage_mat_inl_prev_inc_stored(j_store,:);
end



NewDomSplit(2).history_var_mat_stored=Base_Domain.history_var_mat_stored(j_store,:);
NewDomSplit(2).stress_s1_mat_stored=Base_Domain.stress_s1_mat_stored(j_store,:);
NewDomSplit(2).damage_mat=Base_Domain.damage_mat(j_store,:);
NewDomSplit(2).damage_mat_stored=Base_Domain.damage_mat_stored(j_store,:);
NewDomSplit(2).prop_mat=Base_Domain.prop_mat(j_store,:);
NewDomSplit(2).prop_mat_stored=Base_Domain.prop_mat_stored(j_store,:);


%modify domain 1 to only have properties of the elements remaining
NewDomSplit(1).connect=Base_Domain.connect(:,j_not_store);
NewDomSplit(1).nelnodes=Base_Domain.nelnodes(j_not_store);
NewDomSplit(1).local_strain_mat_stored=Base_Domain.local_strain_mat_stored(j_not_store,:);
NewDomSplit(1).nonlocal_strain_mat_stored= Base_Domain.nonlocal_strain_mat_stored(j_not_store,:);
NewDomSplit(1).history_var_mat=Base_Domain.history_var_mat(j_not_store,:);
NewDomSplit(1).history_var_mat_previousinc= Base_Domain.history_var_mat_previousinc(j_not_store,:);


if SolverID == 3
    NewDomSplit(1).damage_mat_inl=Base_Domain.damage_mat_inl(j_not_store,:);
    NewDomSplit(1).damage_mat_inl_prev_inc=Base_Domain.damage_mat_inl_prev_inc(j_not_store,:);
    NewDomSplit(1).damage_mat_inl_prev_inc_stored=Base_Domain.damage_mat_inl_prev_inc_stored(j_not_store,:);
end


NewDomSplit(1).history_var_mat_stored=Base_Domain.history_var_mat_stored(j_not_store,:);
NewDomSplit(1).stress_s1_mat_stored= Base_Domain.stress_s1_mat_stored(j_not_store,:);
NewDomSplit(1).damage_mat= Base_Domain.damage_mat(j_not_store,:);
NewDomSplit(1).damage_mat_stored= Base_Domain.damage_mat_stored(j_not_store,:);
NewDomSplit(1).prop_mat=Base_Domain.prop_mat(j_not_store,:);
NewDomSplit(1).prop_mat_stored=Base_Domain.prop_mat_stored(j_not_store,:);

l_j_store=length(j_store);
l_j_not_store=length(j_not_store);

%get the number of elements in domain 1 and domain 2
NewDomSplit(2).nelem=l_j_store;
NewDomSplit(1).nelem=l_j_not_store;
% NewDomSplit(1).elident_vec = zeros(NewDomSplit(1).nelem,1);
% NewDomSplit(2).elident_vec = zeros(NewDomSplit(2).nelem,1);
NewDomSplit(1).elident_vec = sparse(NewDomSplit(1).nelem,1);
NewDomSplit(2).elident_vec = sparse(NewDomSplit(2).nelem,1);

%reparametrize the element numbers in bnoth domains
for i =1:l_j_not_store
    NewDomSplit(1).elident_vec(i)=i;
end

for i =1:l_j_store
    NewDomSplit(2).elident_vec(i)=i;
end

%get the coordinate locations of nodes in both domains
coords_loc2=unique(reshape(NewDomSplit(2).connect,[numel(NewDomSplit(2).connect),1]));
NewDomSplit(2).nnodes=length(coords_loc2);
NewDomSplit(2).coords=Base_Domain.coords(:,coords_loc2);

coords_loc1=unique(reshape(NewDomSplit(1).connect,[numel(NewDomSplit(1).connect),1]));
NewDomSplit(1).nnodes=length(coords_loc1);
NewDomSplit(1).coords=Base_Domain.coords(:,coords_loc1);


for i=1:size(NewDomSplit(1).connect,1)
    for j=1:size(NewDomSplit(1).connect,2)
        %go through the connect matrix and and locate their locations
        NewDomSplit(1).connect(i,j)=find(coords_loc1==NewDomSplit(1).connect(i,j));
    end
end

for i=1:size(NewDomSplit(2).connect,1)
    for j=1:size(NewDomSplit(2).connect,2)
        %go through the connect matrix and and locate their locations
        NewDomSplit(2).connect(i,j)=find(coords_loc2==NewDomSplit(2).connect(i,j));
    end
end


%Assumption is that #1 is main domain and #2 is flying domain
NewDomSplit(1).fixnodes=Base_Domain.fixnodes;
NewDomSplit(1).nfix=Base_Domain.nfix;
for i=1:size(NewDomSplit(1).fixnodes,2)
    NewDomSplit(1).fixnodes(1,i)=find(coords_loc1==NewDomSplit(1).fixnodes(1,i));
end
NewDomSplit(2).nfix=0;



%get the number of elements at each node in both domains
NewDomSplit(2).numelem4node=func_Numelem4node(NewDomSplit(2).nnodes, NewDomSplit(2).connect);
NewDomSplit(1).numelem4node=func_Numelem4node(NewDomSplit(1).nnodes, NewDomSplit(1).connect);


% a_vec2 = zeros(2*NewDomSplit(2).nnodes,1);
a_vec2 = sparse(2*NewDomSplit(2).nnodes,1);
for i =1:length(coords_loc2)
    a_vec2(2*i-1)=2*coords_loc2(i)-1;
    a_vec2(2*i)=a_vec2(2*i-1)+1;
end

NewDomSplit(2).Res_F        = Base_Domain.Res_F(a_vec2);
NewDomSplit(2).f_external   = Base_Domain.f_external(a_vec2);
NewDomSplit(2).K            = Base_Domain.K(a_vec2,a_vec2);
NewDomSplit(2).DRdu         = []; % Base_Domain.DRdu(a_vec2,a_vec2);
NewDomSplit(2).dofs         = Base_Domain.dofs(a_vec2);
NewDomSplit(2).dofs_stored  = Base_Domain.dofs_stored(a_vec2);

% a_vec1 = zeros(2*NewDomSplit(1).nnodes,1);
a_vec1 = sparse(2*NewDomSplit(1).nnodes,1);

for i =1:length(coords_loc1)
    a_vec1(2*i-1)=2*coords_loc1(i)-1;
    a_vec1(2*i)=a_vec1(2*i-1)+1;
end

NewDomSplit(1).Res_F=Base_Domain.Res_F(a_vec1);
NewDomSplit(1).f_external=Base_Domain.f_external(a_vec1);
NewDomSplit(1).K=Base_Domain.K(a_vec1,a_vec1);
NewDomSplit(1).DRdu= []; % Base_Domain.DRdu(a_vec1,a_vec1);
NewDomSplit(1).dofs=Base_Domain.dofs(a_vec1);
NewDomSplit(1).dofs_stored=Base_Domain.dofs_stored(a_vec1);


%NewDomSplit(1).boundary_nodes=[];
for i=1:NewDomSplit(1).nnodes
    if(NewDomSplit(1).numelem4node(i)<6)
        NewDomSplit(1).boundary_nodes=[NewDomSplit(1).boundary_nodes,i];
    end
end

for i=1:NewDomSplit(2).nnodes
    if(NewDomSplit(2).numelem4node(i)<4) %this really depends on what mesh you're using. if you have a highly unstructured mesh, it is safe to use 5
        NewDomSplit(2).boundary_nodes=[NewDomSplit(2).boundary_nodes,i];
    end
end

if SolverID == 3

    [NewDomSplit(2).n_hood, NewDomSplit(2).weights] = func_nhood_gausspts(lc, NewDomSplit(2));

end

ndomains=ndomains+1;

[K_pos, fixnodes]= func_bound_pos(NewDomSplit,ndomains);


end