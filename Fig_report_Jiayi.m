figure
f1=plot(H(:,1:11),'r');
hold on
f2=plot(H(:,22:31),'b');
hold on
f3=plot(H(:,42:50),'g');
xlabel('band','fontsize',20);
ylabel('scale','fontsize',20);
legend([f1(1) f2(1) f3(1)],{'Concrete','Asphalt','Colored Structures'});
set(gcf,'outerposition',get(0,'screensize'))
saveas(gcf,'f_3m_spec','epsc');

%% All spectra
figure
f1=plot(H(:,1:11),'r');
xlabel('band','fontsize',20);
ylabel('scale','fontsize',20);
legend(f1(1),'Concrete');
set(gcf,'outerposition',get(0,'screensize'))
saveas(gcf,'f_1con_spec','epsc');

figure
f6=plot(H(:,12:21),'m');
xlabel('band','fontsize',20);
ylabel('scale','fontsize',20);
legend(f6(1),'Red roofs');
set(gcf,'outerposition',get(0,'screensize'))
saveas(gcf,'f_2rr_spec','epsc');

figure
f2=plot(H(:,22:31),'b');
xlabel('band','fontsize',20);
ylabel('scale','fontsize',20);
legend(f2(1),'Asphalt');
set(gcf,'outerposition',get(0,'screensize'))
saveas(gcf,'f_3as_spec','epsc');

figure
f7=plot(H(:,32:41),'y');
xlabel('band','fontsize',20);
ylabel('scale','fontsize',20);
legend(f7(1),'Vegetation');
set(gcf,'outerposition',get(0,'screensize'))
saveas(gcf,'f_4veg_spec','epsc');

figure
f3=plot(H(:,42:50),'g');
xlabel('band','fontsize',20);
ylabel('scale','fontsize',20);
legend(f3(1),'Colored Structures');
set(gcf,'outerposition',get(0,'screensize'))
saveas(gcf,'f_5cs_spec','epsc');


%% plot figure

figure, imshow(rescale(data(:,:,[57,30,20]),1),[])


%% l0g_abun
figure
for p = 1:P   
    subaxis(1,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_g_im(:,:,p,1),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('K_g = 1','fontname','times','fontsize',20)
    end
    colormap jet
    %figname=['f_l0g_abun_Kc1_',num2str(p)];
    set(gcf,'outerposition',get(0,'screensize'))
end
saveas(gcf,'f_l0g_abun_Kc1','epsc');

figure
for p = 1:P
    subaxis(1,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_g_im(:,:,p,2),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('K_g = 2','fontname','times','fontsize',20)
    end
    colormap jet
    set(gcf,'outerposition',get(0,'screensize'))
end
saveas(gcf,'f_l0g_abun_Kc2','epsc');

figure
for p = 1:P   
    subaxis(1,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_g_im(:,:,p,3),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('K_g = 3','fontname','times','fontsize',20)
    end
    colormap jet
    set(gcf,'outerposition',get(0,'screensize'))
end
saveas(gcf,'f_l0g_abun_Kc3','epsc');

figure
for p = 1:P   
    subaxis(1,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)    
    imshow(A_g_im(:,:,p,4),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('K_g = 4','fontname','times','fontsize',20)
    end
    colormap jet
    set(gcf,'outerposition',get(0,'screensize'))
end
saveas(gcf,'f_l0g_abun_Kc4','epsc');

figure
for p = 1:P   
    subaxis(1,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)    
    imshow(A_g_im(:,:,p,5),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('K_g = 5','fontname','times','fontsize',20)
        xlabel('Concrete','fontname','times','fontsize',16)
    elseif p == 2
        xlabel('Red Roofs','fontname','times','fontsize',16)
    elseif p == 3
        xlabel('Asphalt','fontname','times','fontsize',16)
    elseif p == 4
        xlabel('Vegetation','fontname','times','fontsize',16)
    else
        xlabel('Colored Structures','fontname','times','fontsize',16)
    end 
    colormap jet
    set(gcf,'outerposition',get(0,'screensize'))
end
saveas(gcf,'f_l0g_abun_Kc5','epsc');

%set(gcf,'color', 'white')


%% l0gs_abun
% Abundance maps according to groups (MGS-l0)
figure
for p = 1:P    
    subaxis(1,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_mgs_im(:,:,p,1),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        %ylabel(['\tau = ',num2str(tau(1))],'fontname','times','fontsize',20)
        ylabel(['\tau = 1/5'],'fontname','times','fontsize',20)
    end
    colormap jet
    set(gcf,'outerposition',get(0,'screensize'))
end
saveas(gcf,'f_l0gs_abun_tau1','epsc');

figure
for p = 1:P      
    subaxis(1,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_mgs_im(:,:,p,2),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        %ylabel(['\tau = ',num2str(tau(2))],'fontname','times','fontsize',20)
        ylabel(['\tau = 1/10'],'fontname','times','fontsize',20)
    end
    colormap jet
    set(gcf,'outerposition',get(0,'screensize'))
end
saveas(gcf,'f_l0gs_abun_tau2','epsc');

figure
for p = 1:P       
    subaxis(1,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_mgs_im(:,:,p,3),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        %ylabel(['\tau = ',num2str(tau(3))],'fontname','times','fontsize',20)
        ylabel(['\tau = 1/15'],'fontname','times','fontsize',20)
    end
    colormap jet
    set(gcf,'outerposition',get(0,'screensize'))
end
saveas(gcf,'f_l0gs_abun_tau3','epsc');

figure
for p = 1:P       
    subaxis(1,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_mgs_im(:,:,p,4),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        %ylabel(['\tau = ',num2str(tau(4))],'fontname','times','fontsize',20)
        ylabel(['\tau = 1/20'],'fontname','times','fontsize',20)
        xlabel('Concrete','fontname','times','fontsize',16)
    elseif p == 2
        xlabel('Red Roofs','fontname','times','fontsize',16)
    elseif p == 3
        xlabel('Asphalt','fontname','times','fontsize',16)
    elseif p == 4
        xlabel('Vegetation','fontname','times','fontsize',16)
    else
        xlabel('Colored Structures','fontname','times','fontsize',16)
    end 
    colormap jet
    set(gcf,'outerposition',get(0,'screensize'))
end
saveas(gcf,'f_l0gs_abun_tau4','epsc');

%% l0e_abun
figure
for p = 1:P    
    subaxis(1,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_e_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        %ylabel(['\tau = ',num2str(tau(1))],'fontname','times','fontsize',20)
        ylabel(['K_e = 1'],'fontname','times','fontsize',20)
        xlabel('Concrete','fontname','times','fontsize',16)
    elseif p == 2
        xlabel('Red Roofs','fontname','times','fontsize',16)
    elseif p == 3
        xlabel('Asphalt','fontname','times','fontsize',16)
    elseif p == 4
        xlabel('Vegetation','fontname','times','fontsize',16)
    else
        xlabel('Colored Structures','fontname','times','fontsize',16)
    end
    colormap jet
    set(gcf,'outerposition',get(0,'screensize'))
end
saveas(gcf,'f_l0e_abun_Kb1','epsc');

%% abun
% Abundance maps according to groups
figure
for p = 1:P
    subaxis(1,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_g_im(:,:,p,2),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Group (L0)','fontname','times','fontsize',20)
    end
    colormap jet
    set(gcf,'outerposition',get(0,'screensize'))
end
saveas(gcf,'f_l0g_abun','epsc');

figure
for p = 1:P
    subaxis(1,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_e_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Elitist (L0)','fontname','times','fontsize',20)
    end
    colormap jet
    set(gcf,'outerposition',get(0,'screensize'))
end
saveas(gcf,'f_l0e_abun','epsc');

figure
for p = 1:P
    subaxis(1,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_mgs_im(:,:,p,1),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Global (L0)','fontname','times','fontsize',20)
                xlabel('Concrete','fontname','times','fontsize',16)
    elseif p == 2
        xlabel('Red Roofs','fontname','times','fontsize',16)
    elseif p == 3
        xlabel('Asphalt','fontname','times','fontsize',16)
    elseif p == 4
        xlabel('Vegetation','fontname','times','fontsize',16)
    else
        xlabel('Colored Structures','fontname','times','fontsize',16)
    end
    colormap jet
   set(gcf,'outerposition',get(0,'screensize'))
end
saveas(gcf,'f_l0gs_abun','epsc'); 
    
figure
for p = 1:P   
    subaxis(1,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_FCLSU_bundle_corrected_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('FCLSU','fontname','times','fontsize',20)
    end
    colormap jet
    set(gcf,'outerposition',get(0,'screensize'))
end
saveas(gcf,'f_FCLSU_abun','epsc');

figure
for p = 1:P     
    subaxis(1,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_group_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Group','fontname','times','fontsize',20)
    end
    colormap jet
    set(gcf,'outerposition',get(0,'screensize'))
end
saveas(gcf,'f_g_abun','epsc');

figure
for p = 1:P   
    subaxis(1,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_elitist_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Elitist','fontname','times','fontsize',20)
    end
    colormap jet
    set(gcf,'outerposition',get(0,'screensize'))
end
saveas(gcf,'f_e_abun','epsc');

figure
for p = 1:P      
    subaxis(1,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_fractional_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Fractional','fontname','times','fontsize',20)
    end 
    colormap jet
    set(gcf,'outerposition',get(0,'screensize'))
end
saveas(gcf,'f_f_abun','epsc');

%% Qualitative results
% RMSE (data)
figure
imshow(RMSE_g_im(:,:,2),[],'colormap',jet);
colorbar
set(gcf,'outerposition',get(0,'screensize'))
saveas(gcf,'f_RMSE_l0g','epsc');

figure
imshow(RMSE_e_im,[],'colormap',jet);
colorbar
set(gcf,'outerposition',get(0,'screensize'))
saveas(gcf,'f_RMSE_l0e','epsc');

figure
imshow(RMSE_mgs_im(:,:,1),[],'colormap',jet);
colorbar
set(gcf,'outerposition',get(0,'screensize'))
saveas(gcf,'f_RMSE_l0gs','epsc');

figure
imshow(RMSE_FCLSU_bundle_im,[],'colormap', jet)
colorbar
set(gcf,'outerposition',get(0,'screensize'))
saveas(gcf,'f_RMSE_FCLSU','epsc');

figure
imshow(RMSE_group_im,[],'colormap', jet)
colorbar
set(gcf,'outerposition',get(0,'screensize'))
saveas(gcf,'f_RMSE_g','epsc');

figure
imshow(RMSE_elitist_im,[],'colormap', jet)
colorbar
set(gcf,'outerposition',get(0,'screensize'))
saveas(gcf,'f_RMSE_e','epsc');

figure
imshow(RMSE_fractional_im,[],'colormap', jet)
colorbar
set(gcf,'outerposition',get(0,'screensize'))
saveas(gcf,'f_RMSE_f','epsc');


%% [X,Y]£º[8,100]/[60,50] - P1=data(100,8,:) - P2=data(50,60,:)

figure
f1=plot(H(:,1:11),'r');
hold on
f2=plot(H(:,22:31),'b');
hold on
f3=plot(H(:,42:50),'g');
hold on 
f4=plot(Y(:,105*7+100),'k','LineWidth',1.5);
hold on 
f5=plot(Y(:,105*59+50),'c','LineWidth',1.5);
xlabel('band','fontsize',20);
ylabel('scale','fontsize',20);
legend([f1(1) f2(1) f3(1) f4 f5],{'Concrete','Asphalt','Colored Structures','Pixel 1 (8,100)','Pixel 2 (60,50)'});
set(gcf,'outerposition',get(0,'screensize'))
saveas(gcf,'f_spec_P1+P2','epsc');

%% P1-fractional
figure
plot(Y(:,105*7+100),'k','LineWidth',1.5);
hold on
plot(H(:,1:11)*A_fractional_sorted(1:11,105*7+100),'r');
hold on
plot(H(:,12:21)*A_fractional_sorted(12:21,105*7+100),'m');
hold on
plot(H(:,22:31)*A_fractional_sorted(22:31,105*7+100),'b');
hold on
plot(H(:,32:41)*A_fractional_sorted(32:41,105*7+100),'y');
hold on
plot(H(:,42:50)*A_fractional_sorted(42:50,105*7+100),'g');
xlabel('band','fontsize',20);
ylabel('scale','fontsize',20);
legend('Pixel 1 (8,100)','Concrete','Red roofs','Asphalt','Vegetation','Colored Structures');
set(gcf,'outerposition',get(0,'screensize'))
saveas(gcf,'f_spec_unmixP1_f','epsc');

%% P1-global sparsity
figure
plot(Y(:,105*7+100),'k','LineWidth',1.5);
hold on
plot(H(:,1:11)*reshape(A_mgs(1:11,105*7+100,1),11,1),'r');
hold on
plot(H(:,12:21)*reshape(A_mgs(12:21,105*7+100,1),10,1),'m');
hold on
plot(H(:,22:31)*reshape(A_mgs(22:31,105*7+100,1),10,1),'b');
hold on
plot(H(:,32:41)*reshape(A_mgs(32:41,105*7+100,1),10,1),'y');
hold on
plot(H(:,42:50)*reshape(A_mgs(42:50,105*7+100,1),9,1),'g');
xlabel('band','fontsize',20);
ylabel('scale','fontsize',20);
legend('Pixel 1 (8,100)','Concrete','Red roofs','Asphalt','Vegetation','Colored Structures');
set(gcf,'outerposition',get(0,'screensize'))
saveas(gcf,'f_spec_unmixP1_l0gs','epsc');

%% P1-group
figure
plot(Y(:,105*7+100),'k','LineWidth',1.5);
hold on
plot(H(:,1:11)*A_group_sorted(1:11,105*7+100),'r');
hold on
plot(H(:,12:21)*A_group_sorted(12:21,105*7+100),'m');
hold on
plot(H(:,22:31)*A_group_sorted(22:31,105*7+100),'b');
hold on
plot(H(:,32:41)*A_group_sorted(32:41,105*7+100),'y');
hold on
plot(H(:,42:50)*A_group_sorted(42:50,105*7+100),'g');
xlabel('band','fontsize',20);
ylabel('scale','fontsize',20);
legend('Pixel 1 (8,100)','Concrete','Red roofs','Asphalt','Vegetation','Colored Structures');
set(gcf,'outerposition',get(0,'screensize'))
saveas(gcf,'f_spec_unmixP1_g','epsc');

%% P1-group (L0)
figure
plot(Y(:,105*7+100),'k','LineWidth',1.5);
hold on
plot(H(:,1:11)*reshape(A_g(1:11,105*7+100,1),11,1),'r');
hold on
plot(H(:,12:21)*reshape(A_g(12:21,105*7+100,1),10,1),'m');
hold on
plot(H(:,22:31)*reshape(A_g(22:31,105*7+100,1),10,1),'b');
hold on
plot(H(:,32:41)*reshape(A_g(32:41,105*7+100,1),10,1),'y');
hold on
plot(H(:,42:50)*reshape(A_g(42:50,105*7+100,1),9,1),'g');
xlabel('band','fontsize',20);
ylabel('scale','fontsize',20);
legend('Pixel 1 (8,100)','Concrete','Red roofs','Asphalt','Vegetation','Colored Structures');
set(gcf,'outerposition',get(0,'screensize'))
saveas(gcf,'f_spec_unmixP1_l0g','epsc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% P2-fractional
figure
plot(Y(:,105*59+50),'c','LineWidth',1.5);
hold on
plot(H(:,1:11)*A_fractional_sorted(1:11,105*59+50),'r');
hold on
plot(H(:,12:21)*A_fractional_sorted(12:21,105*59+50),'m');
hold on
plot(H(:,22:31)*A_fractional_sorted(22:31,105*59+50),'b');
hold on
plot(H(:,32:41)*A_fractional_sorted(32:41,105*59+50),'y');
hold on
plot(H(:,42:50)*A_fractional_sorted(42:50,105*59+50),'g');
xlabel('band','fontsize',20);
ylabel('scale','fontsize',20);
legend('Pixel 2 (60,50)','Concrete','Red roofs','Asphalt','Vegetation','Colored Structures');
set(gcf,'outerposition',get(0,'screensize'))
saveas(gcf,'f_spec_unmixP2_f','epsc');

%% P2-global sparsity
figure
plot(Y(:,105*59+50),'c','LineWidth',1.5);
hold on
plot(H(:,1:11)*reshape(A_mgs(1:11,105*59+50,1),11,1),'r');
hold on
plot(H(:,12:21)*reshape(A_mgs(12:21,105*59+50,1),10,1),'m');
hold on
plot(H(:,22:31)*reshape(A_mgs(22:31,105*59+50,1),10,1),'b');
hold on
plot(H(:,32:41)*reshape(A_mgs(32:41,105*59+50,1),10,1),'y');
hold on
plot(H(:,42:50)*reshape(A_mgs(42:50,105*59+50,1),9,1),'g');
xlabel('band','fontsize',20);
ylabel('scale','fontsize',20);
legend('Pixel 2 (60,50)','Concrete','Red roofs','Asphalt','Vegetation','Colored Structures');
set(gcf,'outerposition',get(0,'screensize'))
saveas(gcf,'f_spec_unmixP2_l0gs','epsc');

%% P2-group
figure
plot(Y(:,105*59+50),'c','LineWidth',1.5);
hold on
plot(H(:,1:11)*A_group_sorted(1:11,105*59+50),'r');
hold on
plot(H(:,12:21)*A_group_sorted(12:21,105*59+50),'m');
hold on
plot(H(:,22:31)*A_group_sorted(22:31,105*59+50),'b');
hold on
plot(H(:,32:41)*A_group_sorted(32:41,105*59+50),'y');
hold on
plot(H(:,42:50)*A_group_sorted(42:50,105*59+50),'g');
xlabel('band','fontsize',20);
ylabel('scale','fontsize',20);
legend('Pixel 2 (60,50)','Concrete','Red roofs','Asphalt','Vegetation','Colored Structures');
set(gcf,'outerposition',get(0,'screensize'))
saveas(gcf,'f_spec_unmixP2_g','epsc');

%% P2-group (L0)
figure
plot(Y(:,105*59+50),'c','LineWidth',1.5);
hold on
plot(H(:,1:11)*reshape(A_g(1:11,105*59+50,1),11,1),'r');
hold on
plot(H(:,12:21)*reshape(A_g(12:21,105*59+50,1),10,1),'m');
hold on
plot(H(:,22:31)*reshape(A_g(22:31,105*59+50,1),10,1),'b');
hold on
plot(H(:,32:41)*reshape(A_g(32:41,105*59+50,1),10,1),'y');
hold on
plot(H(:,42:50)*reshape(A_g(42:50,105*59+50,1),9,1),'g');
xlabel('band','fontsize',20);
ylabel('scale','fontsize',20);
legend('Pixel 2 (60,50)','Concrete','Red roofs','Asphalt','Vegetation','Colored Structures');
set(gcf,'outerposition',get(0,'screensize'))
saveas(gcf,'f_spec_unmixP2_l0g','epsc');







