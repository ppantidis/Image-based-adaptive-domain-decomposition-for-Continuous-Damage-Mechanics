function [SolverID,TangentID,ncoord,ndof,lc, ...
    increment,inc_success_counter,min_iter,max_iter,max_accept_iter, ...
    loadfactor,dlfactor,dlfactor_incr_threshold,increment_plot_threshold,loadfactor_plot_threshold, ...
    flaglf,countflaglf,incrflag,flagplot, ...
    ndomains,Domain_Vector] = func_read_input_file(infile,model_name)

% ================== READ THE INPUT FILE WITH THE PARAMETERS ==============

% Reads the input text file and stores the variables accordingly
model_parameters_cellarray = textscan(infile,'%s');

cellnumber = 4;

% -------------------------------------------------------------------------
% Solver and Tangent IDs
SolverID    = str2double(model_parameters_cellarray{1}{cellnumber});
cellnumber  = cellnumber + 2;
TangentID   = str2double(model_parameters_cellarray{1}{cellnumber});

% -------------------------------------------------------------------------
% Number of coordinates and degrees of freedom
cellnumber  = cellnumber + 4;
ncoord      = str2double(model_parameters_cellarray{1}{cellnumber});
cellnumber  = cellnumber + 2;
ndof        = str2double(model_parameters_cellarray{1}{cellnumber});

% -------------------------------------------------------------------------
% Characteristic length
cellnumber  = cellnumber + 4;
lc          = str2double(model_parameters_cellarray{1}{cellnumber});

% -------------------------------------------------------------------------
% Adaptive loading parameters
cellnumber                  = cellnumber + 4;
increment                   = str2double(model_parameters_cellarray{1}{cellnumber});
cellnumber                  = cellnumber + 2;
inc_success_counter         = str2double(model_parameters_cellarray{1}{cellnumber});
cellnumber                  = cellnumber + 2;
min_iter                    = str2double(model_parameters_cellarray{1}{cellnumber});
cellnumber                  = cellnumber + 2;
max_iter                    = str2double(model_parameters_cellarray{1}{cellnumber});
cellnumber                  = cellnumber + 2;
max_accept_iter             = str2double(model_parameters_cellarray{1}{cellnumber});
cellnumber                  = cellnumber + 2;
loadfactor                  = str2double(model_parameters_cellarray{1}{cellnumber});
cellnumber                  = cellnumber + 2;
dlfactor                    = str2double(model_parameters_cellarray{1}{cellnumber});
cellnumber                  = cellnumber + 2;
dlfactor_incr_threshold     = str2double(model_parameters_cellarray{1}{cellnumber});
cellnumber                  = cellnumber + 2;
increment_plot_threshold    = str2double(model_parameters_cellarray{1}{cellnumber});
cellnumber                  = cellnumber + 2;
loadfactor_plot_threshold   = str2double(model_parameters_cellarray{1}{cellnumber});
cellnumber                  = cellnumber + 2;
flaglf                      = model_parameters_cellarray{1}{cellnumber};
cellnumber                  = cellnumber + 2;
countflaglf                 = str2double(model_parameters_cellarray{1}{cellnumber});
cellnumber                  = cellnumber + 2;
incrflag                    = str2double(model_parameters_cellarray{1}{cellnumber});
cellnumber                  = cellnumber + 2;
flagplot                    = str2double(model_parameters_cellarray{1}{cellnumber});

% -------------------------------------------------------------------------
% Number of domains
cellnumber  = cellnumber + 4;
ndomains    = str2double(model_parameters_cellarray{1}{cellnumber});

%create Domains
Domain_Vector = domain.empty(ndomains,0);

for i = 1:ndomains
    cellnumber                      = cellnumber + 8;
    Domain_Vector(i)                = domain();
    Domain_Vector(i).ncoord         = ncoord;
    Domain_Vector(i).ndof           = ndof;
    Domain_Vector(i).DomainDamage   = str2double(model_parameters_cellarray{1}{cellnumber});
    cellnumber                      = cellnumber + 3;
    Domain_Vector(i).nprops         = str2double(model_parameters_cellarray{1}{cellnumber});
    Domain_Vector(i).materialprops  = zeros(Domain_Vector(i).nprops,1);

    for j = 1:Domain_Vector(i).nprops
        cellnumber = cellnumber + 2;
        Domain_Vector(i).materialprops(j) = str2double(model_parameters_cellarray{1}{cellnumber});
    end

    cellnumber                  = cellnumber + 3;
    Domain_Vector(i).alpha_val  = str2double(model_parameters_cellarray{1}{cellnumber});
    cellnumber                  = cellnumber + 2;
    Domain_Vector(i).beta_val   = str2double(model_parameters_cellarray{1}{cellnumber});
    cellnumber                  = cellnumber + 2;
    Domain_Vector(i).e_delta    = str2double(model_parameters_cellarray{1}{cellnumber});
    cellnumber                  = cellnumber + 2;
    Domain_Vector(i).dmax       = str2double(model_parameters_cellarray{1}{cellnumber});

    input_name          = strcat(model_name, "_mesh");
    mesh_name           = strcat(input_name,".txt");
    Domain_Vector(i)    = Domain_Vector(i).read_mesh(mesh_name);
    
end



