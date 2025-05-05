function [domain] = func_domain_global_stiffness(domain, IsProj)

% This function calculates the global stiffness, residual and history
% matrices of the domains under consideration

% Include global variables
func_include_flags;

% Initialize variables
% [K, J, DRdu, TangentStiffness, Res_F, history_var_mat, F_int, F_ext] = deal([]);

% Calculate the stiffness matrix from the different domains
% Domain 1: Evaluate stiffness matrix [K] and history variable matrix at the start of the increment
numelem4node                = domain.numelem4node;    % 
dofs                        = domain.dofs;            %
nnodes                      = domain.nnodes;          %
nelem                       = domain.nelem;           %
connect                     = domain.connect;         %
nelnodes                    = domain.nelnodes;        %
maxnodes                    = domain.maxnodes;        %
coords                      = domain.coords;          %
Delastic                    = domain.Delastic;        %
n_hood                      = domain.n_hood;          %
weights                     = domain.weights;         %
ndof                        = domain.ndof;            %
ncoord                      = domain.ncoord;          %
alpha_val                   = domain.alpha_val;       %
beta_val                    = domain.beta_val;        %
e_delta                     = domain.e_delta;         %
dmax                        = domain.dmax;            %
history_var_mat_previousinc = domain.history_var_mat_previousinc; %
DomainDamage                = domain.DomainDamage;

[domain.K, domain.DRdu, domain.Res_F, domain.history_var_mat, domain.damage_mat_inl, domain.gausspoints_prop_mat,domain.nodes_prop_mat] = func_globalstiffness_plain(domain,DomainDamage,dofs,ndof,nnodes,maxnodes,nelnodes,coords,ncoord,nelem,connect,Delastic,history_var_mat_previousinc,alpha_val,beta_val,e_delta,dmax,n_hood,weights,numelem4node,IsProj);

end