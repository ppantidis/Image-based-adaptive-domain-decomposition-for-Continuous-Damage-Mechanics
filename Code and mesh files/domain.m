classdef domain
    % domain: this class defines a domain and its properties

    properties
        nprops                  % No. material parameters
        materialprops           % List of material parameters
        ncoord                  % No. spatial coords (2 for 2D, 3 for 3D)
        ndof                    % No. degrees of freedom per node (2 for 2D, 3 for 3D)
        nnodes                  % No. nodes
        coords                  % coords(i,j)         ith coord of jth node, for i=1..ncoord; j=1..nnode
        nelem                   % No. elements
        maxnodes                % Max no. nodes on any one element
        connect                 % connect(i,j)        List of nodes on the jth element
        nelnodes                % No. nodes on the ith element
        elident_vec             % An integer identifier for the ith element.
        nfix                    % Total no. prescribed displacements
        fixnodes                % fixnodes(i,j)       List of prescribed displacements at nodes

        %fixnodes(1,j) Node number
        %fixnodes(2,j) Displacement component number (1, 2 or 3)
        %fixnodes(3,j) Value of the displacement

        boundary_nodes            %coodinates of the boundary

        Delastic                % Elasticity matrix
        damage_mat              % damage matrix
        damage_mat_stored       % converged damage_mat
        prop_mat                % properties matrix for plotting
        prop_mat_stored         % converged prop_mat
        numelem4node            % number of elements at a specific node
        Res_F                   % internal force vector on domain
        f_external
        K                       % global stiffness Matrix for domain
        DRdu

        % Mazar model Parameters
        alpha_val
        beta_val
        e_delta
        dmax
        DomainDamage
        
        % non-local help matrices
        n_hood
        weights
        damage_mat_inl
        damage_mat_inl_prev_inc
        damage_mat_inl_prev_inc_stored

        dofs                        % guess for current displacements
        dofs_stored                 % converged displacements
        local_strain_mat
        local_strain_mat_stored     % Contains the values of local equivalent strain at each Gauss point
        nonlocal_strain_mat
        nonlocal_strain_mat_stored  % Contains the values of nonlocal equivalent strain at each Gauss point
        history_var_mat
        history_var_mat_stored      % Contains the values of the history variable (damage or kappa) at each Gauss point
        stress_s1_mat
        stress_s1_mat_stored        % Contains the values of 1st principal stress at each Gauss point
        gausspoints_prop_mat
        nodes_prop_mat
        history_var_mat_previousinc

    end

    methods
        function obj = read_mesh(obj, file_name)
            % ===================== READ THE INPUT FILE WITH THE MESH =================
            infile               = fopen(file_name,'r');
            model_mesh_cellarray = textscan(infile,'%s');

            % -------------------------------------------------------------------------
            % Number and coordinates of nodes
            cellnumber      = 2;
            obj.nnodes      = str2double(model_mesh_cellarray{1}{cellnumber});
            cellnumber      = cellnumber + 2;
            obj.coords      = zeros(obj.ncoord,obj.nnodes);
            for i = 1:obj.nnodes
                for j = 1:obj.ncoord
                    obj.coords(j,i) = str2double(model_mesh_cellarray{1}{cellnumber});
                    cellnumber = cellnumber + 1;
                end
            end

            % -------------------------------------------------------------------------
            % Number of elements and their connectivity
            cellnumber      = cellnumber + 1;
            obj.nelem       = str2double(model_mesh_cellarray{1}{cellnumber});
            cellnumber      = cellnumber + 2;
            obj.maxnodes    = str2double(model_mesh_cellarray{1}{cellnumber});
            obj.connect     = zeros(obj.maxnodes,obj.nelem);
            obj.nelnodes    = zeros(obj.nelem,1);
            obj.elident_vec = zeros(obj.nelem,1);
            cellnumber      = cellnumber + 3;

            for i = 1:obj.nelem
                cellnumber          = cellnumber + 1;
                obj.elident_vec(i)  = str2double(model_mesh_cellarray{1}{cellnumber});
                cellnumber          = cellnumber + 1;
                obj.nelnodes(i)     = str2double(model_mesh_cellarray{1}{cellnumber});
                for j = 1:obj.nelnodes(i)
                    cellnumber      = cellnumber + 1;
                    obj.connect(j,i) = str2double(model_mesh_cellarray{1}{cellnumber});
                end
            end

            % -------------------------------------------------------------------------
            % Number of nodes with BCs and prescribed BCs
            cellnumber      = cellnumber + 2;
            obj.nfix        = str2double(model_mesh_cellarray{1}{cellnumber});
            cellnumber      = cellnumber + 3;
            obj.fixnodes    = zeros(3,obj.nfix);

            for i = 1:obj.nfix
                cellnumber          = cellnumber + 1;
                obj.fixnodes(1,i)   = str2double(model_mesh_cellarray{1}{cellnumber});
                cellnumber          = cellnumber + 1;
                obj.fixnodes(2,i)   = str2double(model_mesh_cellarray{1}{cellnumber});
                cellnumber          = cellnumber + 1;
                obj.fixnodes(3,i)   = str2double(model_mesh_cellarray{1}{cellnumber});
            end

            fclose(infile);

            obj.numelem4node = func_Numelem4node(obj.nnodes, obj.connect);
            obj.Delastic     = func_Delastic(obj.materialprops);
            obj.dofs         = zeros(obj.ndof*obj.nnodes,1);
            obj.dofs_stored  = obj.dofs;

            obj.damage_mat          = zeros(obj.nelem,4);
            obj.prop_mat            = zeros(obj.nnodes,7); % 7 is specified based on the number of properties I want to collect.
            obj.damage_mat_stored   = obj.damage_mat;
            obj.prop_mat_stored     = obj.prop_mat;

            % Ensure that the fixnodes matrix lists the prescribed displacements in
            % ascending order of nodes and dofs (as in the global system)
            obj.fixnodes = (sortrows(obj.fixnodes'))';

            % Set up the external and internal force vectors
            obj.Res_F = zeros(obj.ndof*obj.nnodes,1);

            %obj.K = zeros(obj.ndof*obj.nnodes,obj.ndof*obj.nnodes);
            obj.local_strain_mat_stored         = zeros(obj.nelem,4);
            obj.nonlocal_strain_mat_stored      = zeros(obj.nelem,4);
            obj.history_var_mat_stored          = zeros(obj.nelem,4);
            obj.stress_s1_mat_stored            = zeros(obj.nelem,4);
            obj.damage_mat_inl                  = zeros(obj.nelem,4);
            obj.damage_mat_inl_prev_inc_stored  = zeros(obj.nelem,4);
            obj.f_external                      = zeros(obj.ndof*obj.nnodes,1);
        end

        function plotcontours(obj,model_name,inc_success_counter,last_iteration,increment, ...
                loadfactor,Res_F_F_norm,Res_u_norm,solver_string,tangent_string,style)
            func_plotmesh2(obj.coords,obj.connect,obj.nelem,obj.nodes_prop_mat(:,1),obj.nelnodes, ...
                style,model_name,inc_success_counter,last_iteration,increment, ...
                loadfactor,Res_F_F_norm,Res_u_norm,solver_string,tangent_string);
        end
    end
end




