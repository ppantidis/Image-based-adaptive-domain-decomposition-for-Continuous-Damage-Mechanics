classdef element
    %ELEMENT This is a class that defines elements
    %   Detailed explanation goes here

    properties
        lmncoord            %element coordinates
        lmndof              %current displacement for element
        ncoord              %No. spatial coords (2 for 2D, 3 for 3D)
        ndof                %No. degrees of freedom per node (2 for 2D, 3 for 3D)
        Delastic            %Elasticity matrix
        nelnodes            % No. nodes on the element
        ident               %identity number of the element
        damage_mat          %element damage matrix
        k                   %element stiffness matrix
        f_internal          %element internal force vector                        
        prop_mat            %element properties matrix
    end

    methods
        function obj = element(domain,lmn)
            obj.lmncoord = zeros(domain.ncoord,domain.maxnodes);
            obj.lmndof = zeros(domain.ndof,domain.maxnodes);

            for a = 1:domain.nelnodes(lmn)
                for i = 1:domain.ncoord
                    obj.lmncoord(i,a) = domain.coords(i,domain.connect(a,lmn));
                end
                for i = 1:domain.ndof
                    obj.lmndof(i,a) = domain.dofs(domain.ndof*(domain.connect(a,lmn)-1)+i);
                end
            end
            obj.nelnodes = domain.nelnodes(lmn);
            obj.ident = domain.elident(lmn);
            obj.ncoord= domain.ncoord;
            obj.ndof= domain.ndof;
            obj.Delastic = domain.Delastic;
            obj.damage_mat=domain.damage_mat(lmn,:);
            obj.k = zeros(obj.ndof*obj.nelnodes,obj.ndof*obj.nelnodes); 
            obj.f_internal=zeros(obj.ndof*obj.nelnodes,1);
            
        end

    end
end