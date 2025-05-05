function [K,DRdu,dofs,Res_F,f_external,nnodes] = func_assemble(Domain_Vector,ndomains,K_pos_mat,Beta)

func_include_flags;

% Create a bigger K, DRdu, dofs, Res_F, f_external and find nnodes

if TangentID == 1
    DRdu        = [];
    sizeK       = length(Domain_Vector(1).K);
    K           = Domain_Vector(1).K;
    dofs        = Domain_Vector(1).dofs;
    Res_F       = Domain_Vector(1).Res_F;
    f_external  = Domain_Vector(1).f_external;
    nnodes      = Domain_Vector(1).nnodes;


    for j = 2:ndomains
        mainDomSize = sizeK;
        old_K       = K;

        subDomSize  = length(Domain_Vector(j).K);
        sizeK       = sizeK + subDomSize;

        % -----------------------------------------------------------------
        % K = zeros(sizeK, sizeK);
        K = sparse(sizeK, sizeK);

        % -----------------------------------------------------------------

        K(1:mainDomSize,1:mainDomSize)                                               = old_K;
        K(mainDomSize+1:mainDomSize+subDomSize,mainDomSize+1:mainDomSize+subDomSize) = Domain_Vector(j).K;
        

        dofs        = [dofs;Domain_Vector(j).dofs];
        Res_F       = [Res_F;Domain_Vector(j).Res_F];
        f_external  = [f_external;Domain_Vector(j).f_external];

        for i = 1:size(K_pos_mat,2)
            a = K_pos_mat(2*(j-1)-1,i);
            b = K_pos_mat(2*(j-1),i);

            if(a~=0 && b~=0)
                adofx = 2*a-1;
                adofy = 2*a;
                bdofx = 2*b-1;
                bdofy = 2*b;

                K(adofx,adofx)                          = K(adofx,adofx)                         + Beta;
                K(bdofx+mainDomSize,bdofx+mainDomSize)  = K(bdofx+mainDomSize,bdofx+mainDomSize) + Beta;
                K(adofx,bdofx+mainDomSize)              = K(adofx,bdofx+mainDomSize) - Beta;
                K(bdofx+mainDomSize,adofx)              = K(bdofx+mainDomSize,adofx) - Beta;


                K(adofy,adofy)                          = K(adofy,adofy)                            + Beta;
                K(bdofy+mainDomSize,bdofy+mainDomSize)  = K(bdofy+mainDomSize,bdofy+mainDomSize)    + Beta;
                K(adofy,bdofy+mainDomSize)              = K(adofy,bdofy+mainDomSize)                - Beta;
                K(bdofy+mainDomSize,adofy)              = K(bdofy+mainDomSize,adofy)                - Beta;
                
            end
        end
        
        nnodes = nnodes + Domain_Vector(j).nnodes;
    end

elseif TangentID == 2

    K = [];
    if(Domain_Vector(1).DomainDamage == 0)
        sizeDRdu    = length(Domain_Vector(1).K);
        DRdu        = Domain_Vector(1).K;
    else
        sizeDRdu    = length(Domain_Vector(1).DRdu);
        DRdu        = Domain_Vector(1).DRdu;

    end

    dofs        = Domain_Vector(1).dofs;
    Res_F       = Domain_Vector(1).Res_F;
    f_external  = Domain_Vector(1).f_external;
    nnodes      = Domain_Vector(1).nnodes;


    for j = 2:ndomains
        mainDomSize = sizeDRdu;
        old_DRdu    = DRdu;
        subDomSize  = length(Domain_Vector(j).DRdu);
        sizeDRdu    = sizeDRdu + subDomSize;

        DRdu = []; % zeros(sizeDRdu, sizeDRdu);

        DRdu(1:mainDomSize,1:mainDomSize) = old_DRdu;
        if(Domain_Vector(j).DomainDamage==0)
             DRdu(mainDomSize+1:mainDomSize+subDomSize,mainDomSize+1:mainDomSize+subDomSize) = Domain_Vector(j).K;
        else
             DRdu(mainDomSize+1:mainDomSize+subDomSize,mainDomSize+1:mainDomSize+subDomSize) = Domain_Vector(j).DRdu;
        end
       
        for i = 1:size(K_pos_mat,2)
            a = K_pos_mat(2*(j-1)-1,i);
            b = K_pos_mat(2*(j-1),i);

            if(a~=0 && b~=0)
                adofx = 2*a-1;
                adofy = 2*a;
                bdofx = 2*b-1;
                bdofy = 2*b;

                DRdu(adofx,adofx)                           = DRdu(adofx,adofx)                         + Beta;
                DRdu(bdofx+mainDomSize,bdofx+mainDomSize)   = DRdu(bdofx+mainDomSize,bdofx+mainDomSize) + Beta;
                DRdu(adofx,bdofx+mainDomSize)               = DRdu(adofx,bdofx+mainDomSize)             - Beta;
                DRdu(bdofx+mainDomSize,adofx)               = DRdu(bdofx+mainDomSize,adofx)             - Beta;


                DRdu(adofy,adofy)                           = DRdu(adofy,adofy)                         + Beta;
                DRdu(bdofy+mainDomSize,bdofy+mainDomSize)   = DRdu(bdofy+mainDomSize,bdofy+mainDomSize) + Beta;
                DRdu(adofy,bdofy+mainDomSize)               = DRdu(adofy,bdofy+mainDomSize)             - Beta;
                DRdu(bdofy+mainDomSize,adofy)               = DRdu(bdofy+mainDomSize,adofy)             - Beta;

            end
        end


        dofs        = [dofs;Domain_Vector(j).dofs];
        Res_F       = [Res_F;Domain_Vector(j).Res_F];
        f_external  = [f_external;Domain_Vector(j).f_external];
        nnodes      = nnodes + Domain_Vector(j).nnodes;
    end
end






end