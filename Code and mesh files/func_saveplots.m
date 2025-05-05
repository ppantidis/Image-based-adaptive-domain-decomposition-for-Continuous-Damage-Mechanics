function func_saveplots(ndomains,model_name,inc_success_counter,last_iteration,increment, ...
    loadfactor,Res_F_F_norm,Res_u_norm,solver_string,tangent_string,Domain_Vector_graph)

%save stuff====================================
% Initialize figure
% f = figure("Position",[400 200 1100 600]);
% t = tiledlayout(1,2,'TileSpacing','Compact');
% title(t, strcat(model_name + ": Increment #" + num2str(inc_success_counter-1) + " - Iteration #" + num2str(last_iteration) + " - Loadfactor " + num2str(loadfactor)), 'fontweight','bold');
% 
% %figure;
% figure('visible','off')
% 
% nexttile
% hold on; box on;
% colormap(jet);  %https://www.mathworks.com/help/matlab/ref/colormap.html
% 
% 
% for i=1:ndomains
%    
%     Domain_Vector_graph(i).plotcontours(model_name,inc_success_counter,last_iteration,increment, ...
%                 loadfactor,Res_F_F_norm,Res_u_norm,solver_string,tangent_string)
% end
% %axis off
% caxis([0,1])
% 
% %colorbar
% 
% %name2 = strcat('damage_loadfactor_ ', num2str(loadfactor),'.jpg');
% %G=frame2im(getframe(gcf));
% 
% %imwrite(G, name2)
% %exportgraphics(gcf,name2,'Resolution',1000)
% %print(gcf, name2, '-djpeg', '-r1000')
% %saveas(G,name2);
% 
% colorbar
% title("Damage contours", 'fontweight','bold');
% axis equal
% hold off
% 
% % -------------------------------------------------------------------------
% % -------------------------- PLOTTING RESIDUALS ---------------------------
% % -------------------------------------------------------------------------
% nexttile
% semilogy(0:last_iteration,Res_F_F_norm,":ro",'LineWidth',1.1);
% hold on; box on; grid on;
% semilogy(1:last_iteration,Res_u_norm,"--bs",'MarkerFaceColor','b');
% xlim([0, last_iteration]); 
% set(gca, 'XTick', 1:last_iteration); xlabel("No. iterations"); 
% legend(["Internal stresses Residual", "Displacements Residual"]);
% title("Residuals (norm)", 'fontweight','bold');
% hold off
% 
% % -------------------------------------------------------------------------
% saveas(f, strcat(pwd+"\"+model_name + "_Damagecontour_" + solver_string + "_" + tangent_string + "_inc_" + int2str(increment) + ".png"))
% 
% close all 

load custom_colormaps

figure('visible','off'); axis equal; axis off; caxis([0,1]); hold on; colormap(jet);  %https://www.mathworks.com/help/matlab/ref/colormap.html

for i = 1:ndomains
    Domain_Vector_graph(i).plotcontours(model_name,inc_success_counter,last_iteration,increment,loadfactor,Res_F_F_norm,Res_u_norm,solver_string,tangent_string,'interp')
end

name2 = strcat(pwd + "\" + model_name + "_DamContour_" + solver_string + "_inc_" + int2str(increment) + ".png");
% set(gca,'LooseInset',get(gca,'TightInset'));

% Explore printing options in the future: saveas, imwrite, print
exportgraphics(gcf,name2,'Resolution',300)

hold off; close all;




