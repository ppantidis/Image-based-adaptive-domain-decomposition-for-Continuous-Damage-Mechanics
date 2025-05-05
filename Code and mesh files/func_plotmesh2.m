function func_plotmesh2(coords,connect,nelem,wtp,nelnodes,color,model_name,inc_success_counter,last_iteration,increment,loadfactor,Res_F_F_norm,Res_u_norm,solver_string,tangent_string)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =============================== PLOTTING ================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f2D_4 = [1,2,3,4];
%caxis([0 1])
for lmn = 1:nelem
   for i = 1:nelnodes(lmn)
       x(i,1:2) = coords(1:2,connect(i,lmn));
       col(i,1) = wtp(connect(i,lmn));
       %col(i,1) = 0;
   end
   if (nelnodes(lmn)==4)
       patch('Vertices',x,'Faces',f2D_4,'FaceVertexCData',col,'FaceColor','interp','EdgeColor',color);
   else
       disp("Check your nelnodes")
   end
end


end




