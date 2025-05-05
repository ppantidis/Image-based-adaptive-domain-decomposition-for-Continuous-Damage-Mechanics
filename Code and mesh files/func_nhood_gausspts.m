function [neighbours,weights] = func_nhood_gausspts(char_len, domain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======= MATRIX WITH NEIGHBORING GAUSS POINTS (NONLOCAL INTEGRAL) ========
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include global variables
func_include_flags;

% -------------------------------------------------------------------------
lmncoord = zeros(domain.ncoord,domain.maxnodes);
x_coord  = zeros(domain.nelem*4,2);


for lmn = 1:domain.nelem
    
    % Extract coords of nodes
    for a = 1:domain.nelnodes(lmn)
        for i = 1:domain.ncoord
            lmncoord(i,a) = domain.coords(i,domain.connect(a,lmn));
        end
    end
    
    % This is the same as lmn. I just left it in because I had already
    % written it
    ident = domain.elident_vec(lmn);
    
    % Setting up:
    npoints = func_numberofintegrationpoints;  % Number of integration points
    xilist  = func_integrationpoints;          % Positions of integration points
    
    % Loop over the integration points
    for intpt = 1:npoints
        
        % Compute shape functions && derivatives wrt local coords
        for i = 1:domain.ncoord
            xi(i) = xilist(i,intpt);
        end
        
        N = func_shapefunctions(xi);   % Shape function for each integration point
        
        x_coord(4 * ident + intpt - 4,:) = N' * lmncoord';% Coordinates of integration point
        
    end
    
end

for i = 1:npoints*domain.nelem
    cnt4nopnear = 0; % counter for no points near
    count = 1;
    % If a point is within the radius of the characteristic length add it to
    %the neighbour matrix and calculate its weight
    for j = 1:npoints*domain.nelem
        diff_distance = norm(x_coord(i,:) - x_coord(j,:));
        
        if diff_distance <= char_len
            neighbours(i,count) = j;
            weights(i,count)    = exp(-diff_distance / (char_len * char_len));
            count               = count + 1;
            cnt4nopnear         = 1;
        end
        
    end
    
    % If there are no points near a point, make it a neighbour of itself
    %with a weight of 1
    if cnt4nopnear == 0
        neighbours(i,count) = i;
        weights(i,count)    = 1;
    end
   
    
end

end


