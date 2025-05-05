function [K_pos, fixnodes]= func_bound_pos(Domain_Vector,ndomains)
%This function defines the connections of the boundaries relative to each other

%if there is only one domain, then this matrix is not important, so return
%0
if ndomains ==1
    K_pos=0;
    fixnodes=Domain_Vector(1).fixnodes;
    return

elseif ndomains > 1
    
    %make a matrix which has size ndomains by max # of nodes at a boundary
    %this is just for speed of allocation purposes
    max_nodes_bound=0;

    for i=2:ndomains
        max_nodes_bound=max(max_nodes_bound, length(Domain_Vector(i).boundary_nodes));
    end

    % K_pos = zeros(2*(ndomains-1),max_nodes_bound);
    K_pos = sparse(2*(ndomains-1),max_nodes_bound);    


    %go through every subdomain from #3 to the end and find how the
    %boundaries connect each other
    for num_sub_domain=2:ndomains
        counter=1;
        for i=1:length(Domain_Vector(num_sub_domain).boundary_nodes)
            
            %boundary coords on the current subdomain
            check=Domain_Vector(num_sub_domain).coords(:,Domain_Vector(num_sub_domain).boundary_nodes(i));

            for j=1:length(Domain_Vector(1).boundary_nodes)
                %boundary coords on subdomains which have already been
                %joined
                check_coords2=Domain_Vector(1).coords(:,Domain_Vector(1).boundary_nodes(j));
                if norm(check(1)-check_coords2(1))<=10^-12 && norm(check(2)-check_coords2(2))<=10^-12
                    K_pos(2*(num_sub_domain-1)-1,counter)=Domain_Vector(1).boundary_nodes(j);
                    K_pos(2*(num_sub_domain-1),counter)=Domain_Vector(num_sub_domain).boundary_nodes(i);
                    counter=counter+1;
                    break
                end
            end
        end
        
        %every time you are done comparing two domains, join them
        Domain_Vector(1).coords=cat(2,Domain_Vector(1).coords,Domain_Vector(num_sub_domain).coords);
        Domain_Vector(num_sub_domain).boundary_nodes=length(Domain_Vector(1).nnodes) +...
            Domain_Vector(num_sub_domain).boundary_nodes;
        Domain_Vector(1).boundary_nodes=cat(2,Domain_Vector(1).boundary_nodes,Domain_Vector(num_sub_domain).boundary_nodes);
        
        %combine fixnodes in case both domains have their own fixnodes
        if(~isempty(Domain_Vector(num_sub_domain).fixnodes))
            Domain_Vector(num_sub_domain).fixnodes(1,:)=length(Domain_Vector(1).nnodes)+Domain_Vector(num_sub_domain).fixnodes(1,:);
            Domain_Vector(1).fixnodes=cat(2,Domain_Vector(1).fixnodes,Domain_Vector(num_sub_domain).fixnodes);
        end



    end
    fixnodes=Domain_Vector(1).fixnodes;

end
end