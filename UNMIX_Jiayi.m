clc
clear all
close all
dbstop if error

seed = 1; % super duper important

rng(seed) % super duper important

load('real_data_1.mat');

figure, imshow(rescale(data(:,:,[57,30,20]),1),[])

% data =rescale(x20100628_gipsa_area2(300:499,1:200,:),1) ;

EMs = figure;
imshow(rescale(data(:,:,[57,30,20]),1),[])
hold on

[m,n,L] = size(data);

%% Extract endmember candidates and cluster them into bundles

%  load bundle_test_paper_test.mat

Y = reshape(data,m*n,L)';

imr = Y;

P = 5;
% load bundle_test_paper.mat

bundle_nbr = 10;
percent = 10;
clustering = 'kmeans';
% clustering = 'spectral_clustering';
% clustering = 'nystrom'; % just to test function

% groups = kmeans(imr',P,'distance','cosine');
% bundle = imr;

[groups, bundle] = batchvca_modif(Y, P, bundle_nbr, percent, clustering );

%% Reshape dictionary
dic=sortrows([groups';bundle]',1);
dic=dic';
spec_group=dic(1,:);
H=dic(2:end,:);
nbelements=tabulate(spec_group);
nbelements=nbelements(:,2);

%% Group_L0
type='group';
A_g=zeros(length(spec_group),m*n,P);
t_g=zeros(1,P);
for K=1:P
    [A,est,t] = mip_l0(H,Y,[],K,[],nbelements,spec_group,type);
    A_g(:,:,K)=A;
    t_g(K)=sum(t);
end

% sum the abundances within each class (bundle2global function)
A_g_final=zeros(P,m*n,P);
A_g_im=zeros(m,n,P,P);
A_g_full_im=zeros(m,n,P*bundle_nbr,P);
H_g=zeros(L,m*n,P);
RMSE_g=zeros(P,m*n);
RMSE_g_im=zeros(m,n,P);
for K=1:P
    [A_g_final(:,:,K),S_g] = bundle2global(A_g(:,:,K),H,spec_group');
    A_g_im(:,:,:,K) = reshape(A_g_final(:,:,K)',m,n,P);
    A_g_full_im(:,:,:,K) = reshape(A_g(:,:,K)',m,n,P*bundle_nbr);
    H_g(:,:,K) = H*A_g(:,:,K); 
    RMSE_g(K,:) = sqrt(1/L*sum((H_g(:,:,K)-Y).^2,1));
    RMSE_g_im(:,:,K) = reshape(RMSE_g(K,:),m,n);
end

SAM_g=zeros(m*n,P);
SAM_g_im=zeros(m,n,P);
for K=1:P
    for k = 1:m*n
        SAM_g(k,K) = 180/pi*real(acos((Y(:,k)'*H_g(:,k,K))/(norm(Y(:,k))*norm(H_g(:,k,K)))));
    end
    SAM_g_im(:,:,K)=reshape(SAM_g(:,K),m,n);
end

% save('Jiayi_G.mat')

%% Elitist_L0
type='elitist';
[A_e,est_e,t_e] = mip_l0(H,Y,P,[],[],nbelements,spec_group,type);
t_e=sum(t_e);

[A_e_final,S_e] = bundle2global(A_e,H,spec_group');
H_e = H*A_e; 
RMSE_e = sqrt(1/L*sum((H_e-Y).^2,1));       
RMSE_e_im = reshape(RMSE_e,m,n);

SAM_e = zeros(m*n,1);
for k = 1:m*n
        SAM_e(k) = 180/pi*real(acos((Y(:,k)'*H_e(:,k))/(norm(Y(:,k))*norm(H_e(:,k)))));
end
SAM_e_im=reshape(SAM_e,m,n);

A_e_im = reshape(A_e_final',m,n,P);
A_e_full_im = reshape(A_e',m,n,P*bundle_nbr);

% save('Jiayi_E.mat')

%% Global_Sparsity_L0
Kb=5:5:20;
Kc=2;
K=length(Kb);
type='gs';
A_gs=zeros(length(spec_group),m*n,K);
t_gs=zeros(1,K);
for i=1:K
    [A,est,t] = mip_l0(H,Y,Kb(i),Kc,[],nbelements,spec_group,type);
    A_gs(:,:,i)=A;
    t_gs(i)=sum(t);
end

A_gs_final=zeros(P,m*n,K);
A_gs_im=zeros(m,n,P,K);
A_gs_full_im=zeros(m,n,P*bundle_nbr,K);
H_gs=zeros(L,m*n,K);
RMSE_gs=zeros(K,m*n);
RMSE_gs_im=zeros(m,n,K);
for i=1:K
    [A_gs_final(:,:,i),S_gs] = bundle2global(A_gs(:,:,i),H,spec_group');
    A_gs_im(:,:,:,i) = reshape(A_gs_final(:,:,i)',m,n,P);
    A_gs_full_im(:,:,:,i) = reshape(A_gs(:,:,i)',m,n,P*bundle_nbr);
    H_gs(:,:,i) = H*A_gs(:,:,i); 
    RMSE_gs(i,:) = sqrt(1/L*sum((H_gs(:,:,i)-Y).^2,1));
    RMSE_gs_im(:,:,i) = reshape(RMSE_gs(i,:),m,n);
end

SAM_gs=zeros(m*n,K);
SAM_gs_im=zeros(m,n,K);
for i=1:K
    for k = 1:m*n
        SAM_gs(k,i) = 180/pi*real(acos((Y(:,k)'*H_gs(:,k,i))/(norm(Y(:,k))*norm(H_gs(:,k,i)))));
    end
    SAM_gs_im(:,:,i)=reshape(SAM_gs(:,i),m,n);
end

save(['Jiayi_GS_Kc(',num2str(Kc),').mat'])

% save(['Jiayi_GS_Kb(',num2str(Kb),')_Kc(',num2str(Kc),').mat'])

%% Modified_Global_Sparsity_L0
Kb=5:5:20;
Kc=3;
tau=1./Kb;
K=length(tau);
type='mgs';
A_mgs=zeros(length(spec_group),m*n,K);
t_mgs=zeros(1,K);
for i=1:K
    [A,est,t] = mip_l0(H,Y,[],Kc,tau(i),nbelements,spec_group,type);
    A_mgs(:,:,i)=A;
    t_mgs(i)=sum(t);
end

A_mgs_final=zeros(P,m*n,K);
A_mgs_im=zeros(m,n,P,K);
A_mgs_full_im=zeros(m,n,P*bundle_nbr,K);
H_mgs=zeros(L,m*n,K);
RMSE_mgs=zeros(K,m*n);
RMSE_mgs_im=zeros(m,n,K);
for i=1:K
    [A_mgs_final(:,:,i),S_gs] = bundle2global(A_mgs(:,:,i),H,spec_group');
    A_mgs_im(:,:,:,i) = reshape(A_mgs_final(:,:,i)',m,n,P);
    A_mgs_full_im(:,:,:,i) = reshape(A_mgs(:,:,i)',m,n,P*bundle_nbr);
    H_mgs(:,:,i) = H*A_mgs(:,:,i); 
    RMSE_mgs(i,:) = sqrt(1/L*sum((H_mgs(:,:,i)-Y).^2,1));
    RMSE_mgs_im(:,:,i) = reshape(RMSE_mgs(i,:),m,n);
end

SAM_mgs=zeros(m*n,K);
SAM_mgs_im=zeros(m,n,K);
for i=1:K
    for k = 1:m*n
        SAM_mgs(k,i) = 180/pi*real(acos((Y(:,k)'*H_mgs(:,k,i))/(norm(Y(:,k))*norm(H_mgs(:,k,i)))));
    end
    SAM_mgs_im(:,:,i)=reshape(SAM_mgs(:,i),m,n);
end

save(['Jiayi_MGS_Kc(',num2str(Kc),').mat'])

% save(['Jiayi_MGS_Kc(',num2str(Kc),')_tau(',num2str(tau),').mat'])

%% 
% save('Jiayi_data.mat')
% save('Jiayi_data_Iter1000.mat')

%% Qualitative results - Group(L0-norm)
% RMSE (Group_L0)
figure
suptitle('The average data RMSE - Group(L0)');
subplot(2,3,1)
imshow(RMSE_g_im(:,:,1),[],'colormap',jet);
title(['K=1, mean: ',num2str(mean(RMSE_g(1,:)))]);
colorbar
subplot(2,3,2)
imshow(RMSE_g_im(:,:,2),[],'colormap',jet);
title(['K=2, mean: ',num2str(mean(RMSE_g(2,:)))]);
colorbar
subplot(2,3,3)
imshow(RMSE_g_im(:,:,3),[],'colormap', jet)
title(['K=3, mean: ',num2str(mean(RMSE_g(3,:)))])
colorbar
subplot(2,3,4)
imshow(RMSE_g_im(:,:,4),[],'colormap', jet)
title(['K=4, mean: ',num2str(mean(RMSE_g(4,:)))])
colorbar
subplot(2,3,5)
imshow(RMSE_g_im(:,:,5),[],'colormap', jet)
title(['K=5, mean: ',num2str(mean(RMSE_g(5,:)))])
colorbar
set(gcf,'color', 'white');

%%
% SAM (Group_L0)
figure
suptitle('The average SAM - Group(L0)');
subplot(2,3,1)
imshow(SAM_g_im(:,:,1),[],'colormap',jet);
title(['K=1, mean: ',num2str(mean(SAM_g(:,1)))]);
colorbar
subplot(2,3,2)
imshow(SAM_g_im(:,:,2),[],'colormap',jet);
title(['K=2, mean: ',num2str(mean(SAM_g(:,2)))]);
colorbar
subplot(2,3,3)
imshow(SAM_g_im(:,:,3),[],'colormap', jet)
title(['K=3, mean: ',num2str(mean(SAM_g(:,3)))])
colorbar
subplot(2,3,4)
imshow(SAM_g_im(:,:,4),[],'colormap', jet)
title(['K=4, mean: ',num2str(mean(SAM_g(:,4)))])
colorbar
subplot(2,3,5)
imshow(SAM_g_im(:,:,5),[],'colormap', jet)
title(['K=5, mean: ',num2str(mean(SAM_g(:,5)))])
colorbar
set(gcf,'color', 'white');

%% 
% Abundance maps according to groups (Group-l0)
figure
suptitle('Abundance maps - Group(L0)')
for p = 1:P
    
    subaxis(5,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_g_im(:,:,p,1),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('K=1','fontname','times','fontsize',12)
    end
    colormap jet
    
    subaxis(5,P,P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_g_im(:,:,p,2),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('K=2','fontname','times','fontsize',12)
    end
    colormap jet
    
    subaxis(5,P,2*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_g_im(:,:,p,3),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('K=3','fontname','times','fontsize',12)
    end
    colormap jet
    
    subaxis(5,P,3*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_g_im(:,:,p,4),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('K=4','fontname','times','fontsize',12)
    end
    colormap jet
    
    subaxis(5,P,4*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_g_im(:,:,p,5),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('K=5','fontname','times','fontsize',12)
        xlabel('Concrete','fontname','times','fontsize',12)
    elseif p == 2
        xlabel('Red Roofs','fontname','times','fontsize',12)
    elseif p == 3
        xlabel('Asphalt','fontname','times','fontsize',12)
    elseif p == 4
        xlabel('Vegetation','fontname','times','fontsize',12)
    else
        xlabel('Colored Structures','fontname','times','fontsize',12)
    end 
    colormap jet
end
set(gcf,'color', 'white')

%% 
% Abundance maps according to endmembers (Group-l0)
for  i = 1:P
    EMs_L0_cons = find(spec_group == i);
    
    figure
    suptitle(['Abundance maps of group ',num2str(i),' - Group(L0)'])
    for p = 1:length(EMs_L0_cons)
        
        subaxis(5,length(EMs_L0_cons),p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
        imshow(A_g_full_im(:,:,EMs_L0_cons(p),1),[],'colormap',jet)
        set(gca,'clim',[0,1])
        if p == 1
            ylabel('K=1','fontname','times','fontsize',12)
        end
        
        subaxis(5,length(EMs_L0_cons),length(EMs_L0_cons)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
        imshow(A_g_full_im(:,:,EMs_L0_cons(p),2),[],'colormap',jet)
        set(gca,'clim',[0,1])
        if p == 1
            ylabel('K=2','fontname','times','fontsize',12)
        end
        
        subaxis(5,length(EMs_L0_cons),2*length(EMs_L0_cons)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
        imshow(A_g_full_im(:,:,EMs_L0_cons(p),3),[],'colormap',jet)
        set(gca,'clim',[0,1])
        if p == 1
            ylabel('K=3','fontname','times','fontsize',12)
        end
        
        subaxis(5,length(EMs_L0_cons),3*length(EMs_L0_cons)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
        imshow(A_g_full_im(:,:,EMs_L0_cons(p),4),[],'colormap',jet)
        set(gca,'clim',[0,1])
        if p == 1
            ylabel('K=4','fontname','times','fontsize',12)
        end
        
        subaxis(5,length(EMs_L0_cons),4*length(EMs_L0_cons)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
        imshow(A_g_full_im(:,:,EMs_L0_cons(p),5),[],'colormap',jet)
        set(gca,'clim',[0,1])
        if p == 1
            ylabel('K=5','fontname','times','fontsize',12)
        end
        set(gcf,'colormap', jet)
    end
    set(gcf,'color', 'white')
end


%% Qualitative results - Global Sparsity(L0-norm)
% RMSE (GS_L0)
figure
suptitle('The average data RMSE - Global Sparsity(L0)');
subplot(2,2,1)
imshow(RMSE_gs_im(:,:,1),[],'colormap',jet);
title(['Kb=5, mean: ',num2str(mean(RMSE_gs(1,:)))]);
colorbar
subplot(2,2,2)
imshow(RMSE_gs_im(:,:,2),[],'colormap',jet);
title(['Kb=10, mean: ',num2str(mean(RMSE_gs(2,:)))]);
colorbar
subplot(2,2,3)
imshow(RMSE_gs_im(:,:,3),[],'colormap', jet)
title(['Kb=15, mean: ',num2str(mean(RMSE_gs(3,:)))])
colorbar
subplot(2,2,4)
imshow(RMSE_gs_im(:,:,4),[],'colormap', jet)
title(['Kb=20, mean: ',num2str(mean(RMSE_gs(4,:)))])
colorbar
set(gcf,'color', 'white');

%%
% SAM (GS_L0)
figure
suptitle('The average SAM - Global Sparsity(L0)');
subplot(2,2,1)
imshow(SAM_gs_im(:,:,1),[],'colormap',jet);
title(['Kb=5, mean: ',num2str(mean(SAM_gs(:,1)))]);
colorbar
subplot(2,2,2)
imshow(SAM_gs_im(:,:,2),[],'colormap',jet);
title(['Kb=10, mean: ',num2str(mean(SAM_gs(:,2)))]);
colorbar
subplot(2,2,3)
imshow(SAM_gs_im(:,:,3),[],'colormap', jet)
title(['Kb=15, mean: ',num2str(mean(SAM_gs(:,3)))])
colorbar
subplot(2,2,4)
imshow(SAM_gs_im(:,:,4),[],'colormap', jet)
title(['Kb=20, mean: ',num2str(mean(SAM_gs(:,4)))])
colorbar
set(gcf,'color', 'white');

%% 
% Abundance maps according to groups (GS-l0)
figure
suptitle('Abundance maps - Global Sparsity(L0)')
for p = 1:P
    
    subaxis(4,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_gs_im(:,:,p,1),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Kb=5','fontname','times','fontsize',12)
    end
    colormap jet
    
    subaxis(4,P,P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_gs_im(:,:,p,2),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Kb=10','fontname','times','fontsize',12)
    end
    colormap jet
    
    subaxis(4,P,2*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_gs_im(:,:,p,3),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Kb=15','fontname','times','fontsize',12)
    end
    colormap jet
    
    subaxis(4,P,3*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_gs_im(:,:,p,4),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Kb=20','fontname','times','fontsize',12)
        xlabel('Concrete','fontname','times','fontsize',12)
    elseif p == 2
        xlabel('Red Roofs','fontname','times','fontsize',12)
    elseif p == 3
        xlabel('Asphalt','fontname','times','fontsize',12)
    elseif p == 4
        xlabel('Vegetation','fontname','times','fontsize',12)
    else
        xlabel('Colored Structures','fontname','times','fontsize',12)
    end 
    colormap jet
    
end
set(gcf,'color', 'white')

%% 
% Abundance maps according to endmembers (GS-l0)
for  i = 1:P
    EMs_L0_cons = find(spec_group == i);
    
    figure
    suptitle(['Abundance maps of group ',num2str(i),' - Global Sparsity(L0)'])
    for p = 1:length(EMs_L0_cons)
        
        subaxis(4,length(EMs_L0_cons),p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
        imshow(A_gs_full_im(:,:,EMs_L0_cons(p),1),[],'colormap',jet)
        set(gca,'clim',[0,1])
        if p == 1
            ylabel('Kb=5','fontname','times','fontsize',12)
        end
        
        subaxis(4,length(EMs_L0_cons),length(EMs_L0_cons)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
        imshow(A_gs_full_im(:,:,EMs_L0_cons(p),2),[],'colormap',jet)
        set(gca,'clim',[0,1])
        if p == 1
            ylabel('Kb=10','fontname','times','fontsize',12)
        end
        
        subaxis(4,length(EMs_L0_cons),2*length(EMs_L0_cons)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
        imshow(A_gs_full_im(:,:,EMs_L0_cons(p),3),[],'colormap',jet)
        set(gca,'clim',[0,1])
        if p == 1
            ylabel('Kb=15','fontname','times','fontsize',12)
        end
        
        subaxis(4,length(EMs_L0_cons),3*length(EMs_L0_cons)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
        imshow(A_gs_full_im(:,:,EMs_L0_cons(p),4),[],'colormap',jet)
        set(gca,'clim',[0,1])
        if p == 1
            ylabel('Kb=20','fontname','times','fontsize',12)
        end
        set(gcf,'colormap', jet)
    end
    set(gcf,'color', 'white')
end


%% Qualitative results - Modified Global Sparsity(L0-norm)
% RMSE (MGS_L0)
figure
suptitle('The average data RMSE - MGS(L0)');
subplot(2,2,1)
imshow(RMSE_mgs_im(:,:,1),[],'colormap',jet);
title(['tau=',num2str(tau(1)),' mean: ',num2str(mean(RMSE_mgs(1,:)))]);
colorbar
subplot(2,2,2)
imshow(RMSE_mgs_im(:,:,2),[],'colormap',jet);
title(['tau=',num2str(tau(2)),' mean: ',num2str(mean(RMSE_mgs(2,:)))]);
colorbar
subplot(2,2,3)
imshow(RMSE_mgs_im(:,:,3),[],'colormap', jet)
title(['tau=',num2str(tau(3)),' mean: ',num2str(mean(RMSE_mgs(3,:)))]);
colorbar
subplot(2,2,4)
imshow(RMSE_mgs_im(:,:,4),[],'colormap', jet)
title(['tau=',num2str(tau(4)),' mean: ',num2str(mean(RMSE_mgs(4,:)))]);
colorbar
set(gcf,'color', 'white');

%%
% SAM (MGS_L0)
figure
suptitle('The average SAM - MGS(L0)');
subplot(2,2,1)
imshow(SAM_mgs_im(:,:,1),[],'colormap',jet);
title(['tau=',num2str(tau(1)),' mean: ',num2str(mean(SAM_mgs(:,1)))]);
colorbar
subplot(2,2,2)
imshow(SAM_mgs_im(:,:,2),[],'colormap',jet);
title(['tau=',num2str(tau(2)),' mean: ',num2str(mean(SAM_mgs(:,2)))]);
colorbar
subplot(2,2,3)
imshow(SAM_mgs_im(:,:,3),[],'colormap', jet)
title(['tau=',num2str(tau(3)),' mean: ',num2str(mean(SAM_mgs(:,3)))])
colorbar
subplot(2,2,4)
imshow(SAM_mgs_im(:,:,4),[],'colormap', jet)
title(['tau=',num2str(tau(4)),' mean: ',num2str(mean(SAM_mgs(:,4)))])
colorbar
set(gcf,'color', 'white');

%% 
% Abundance maps according to groups (MGS-l0)
figure
suptitle('Abundance maps - MGS(L0)')
for p = 1:P
    
    subaxis(4,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_mgs_im(:,:,p,1),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel(['tau=',num2str(tau(1))],'fontname','times','fontsize',12)
    end
    colormap jet
    
    subaxis(4,P,P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_mgs_im(:,:,p,2),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel(['tau=',num2str(tau(2))],'fontname','times','fontsize',12)
    end
    colormap jet
    
    subaxis(4,P,2*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_mgs_im(:,:,p,3),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel(['tau=',num2str(tau(3))],'fontname','times','fontsize',12)
    end
    colormap jet
    
    subaxis(4,P,3*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_mgs_im(:,:,p,4),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel(['tau=',num2str(tau(4))],'fontname','times','fontsize',12)
        xlabel('Concrete','fontname','times','fontsize',12)
    elseif p == 2
        xlabel('Red Roofs','fontname','times','fontsize',12)
    elseif p == 3
        xlabel('Asphalt','fontname','times','fontsize',12)
    elseif p == 4
        xlabel('Vegetation','fontname','times','fontsize',12)
    else
        xlabel('Colored Structures','fontname','times','fontsize',12)
    end 
    colormap jet
    
end
set(gcf,'color', 'white')

%% 
% Abundance maps according to endmembers (MGS-l0)
for  i = 1:P
    EMs_L0_cons = find(spec_group == i);
    
    figure
    suptitle(['Abundance maps of group ',num2str(i),' - MGS(L0)'])
    for p = 1:length(EMs_L0_cons)
        
        subaxis(4,length(EMs_L0_cons),p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
        imshow(A_mgs_full_im(:,:,EMs_L0_cons(p),1),[],'colormap',jet)
        set(gca,'clim',[0,1])
        if p == 1
            ylabel(['tau=',num2str(tau(1))],'fontname','times','fontsize',12)
        end
        
        subaxis(4,length(EMs_L0_cons),length(EMs_L0_cons)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
        imshow(A_mgs_full_im(:,:,EMs_L0_cons(p),2),[],'colormap',jet)
        set(gca,'clim',[0,1])
        if p == 1
            ylabel(['tau=',num2str(tau(2))],'fontname','times','fontsize',12)
        end
        
        subaxis(4,length(EMs_L0_cons),2*length(EMs_L0_cons)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
        imshow(A_mgs_full_im(:,:,EMs_L0_cons(p),3),[],'colormap',jet)
        set(gca,'clim',[0,1])
        if p == 1
            ylabel(['tau=',num2str(tau(3))],'fontname','times','fontsize',12)
        end
        
        subaxis(4,length(EMs_L0_cons),3*length(EMs_L0_cons)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
        imshow(A_mgs_full_im(:,:,EMs_L0_cons(p),4),[],'colormap',jet)
        set(gca,'clim',[0,1])
        if p == 1
            ylabel(['tau=',num2str(tau(4))],'fontname','times','fontsize',12)
        end
        set(gcf,'colormap', jet)
    end
    set(gcf,'color', 'white')
end

%% Compare_GS_MGS
% Abundance maps according to groups (GS-MGS-l0)
figure
suptitle('Abundance maps - GS/MGS(L0)')
for p = 1:P
    
    subaxis(8,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_gs_im(:,:,p,1),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Kb=5','fontname','times','fontsize',12)
    end
    colormap jet
    
    subaxis(8,P,2*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_gs_im(:,:,p,2),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Kb=10','fontname','times','fontsize',12)
    end
    colormap jet
    
    subaxis(8,P,4*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_gs_im(:,:,p,3),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Kb=15','fontname','times','fontsize',12)
    end
    colormap jet
    
    subaxis(8,P,6*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_gs_im(:,:,p,4),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Kb=20','fontname','times','fontsize',12)
    end 
    colormap jet
    
    subaxis(8,P,P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_mgs_im(:,:,p,1),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel(['tau=',num2str(tau(1))],'fontname','times','fontsize',12)
    end
    colormap jet
    
    subaxis(8,P,3*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_mgs_im(:,:,p,2),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel(['tau=',num2str(tau(2))],'fontname','times','fontsize',12)
    end
    colormap jet
    
    subaxis(8,P,5*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_mgs_im(:,:,p,3),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel(['tau=',num2str(tau(3))],'fontname','times','fontsize',12)
    end
    colormap jet
    
    subaxis(8,P,7*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_mgs_im(:,:,p,4),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel(['tau=',num2str(tau(4))],'fontname','times','fontsize',12)
        xlabel('Concrete','fontname','times','fontsize',12)
    elseif p == 2
        xlabel('Red Roofs','fontname','times','fontsize',12)
    elseif p == 3
        xlabel('Asphalt','fontname','times','fontsize',12)
    elseif p == 4
        xlabel('Vegetation','fontname','times','fontsize',12)
    else
        xlabel('Colored Structures','fontname','times','fontsize',12)
    end 
    colormap jet
    
end
set(gcf,'color', 'white')


%% Qualitative results
load('Lucas_data.mat');  % save from houston_paper

% RMSE (data)
figure
suptitle('The average data RMSE');
subplot(2,4,1)
imshow(RMSE_g_im(:,:,2),[],'colormap',jet);
title(['RMSE group(L0), mean: ',num2str(mean(RMSE_gs(2,:)))]);
colorbar
subplot(2,4,2)
imshow(RMSE_e_im,[],'colormap',jet);
title(['RMSE elitist(L0), mean: ',num2str(mean(RMSE_e(:)))]);
colorbar
subplot(2,4,3)
imshow(RMSE_mgs_im(:,:,1),[],'colormap',jet);
title(['RMSE global(L0), mean: ',num2str(mean(RMSE_mgs(1,:)))]);
colorbar
subplot(2,4,4)
imshow(RMSE_FCLSU_bundle_im,[],'colormap', jet)
title(['RMSE FCLSU bundle, mean: ',num2str(mean(RMSE_FCLSU_bundle(:)))])
colorbar
subplot(2,4,5)
imshow(RMSE_group_im,[],'colormap', jet)
title(['RMSE group, mean: ',num2str(mean(RMSE_group(:)))])
colorbar
subplot(2,4,6)
imshow(RMSE_elitist_im,[],'colormap', jet)
title(['RMSE elitist, mean: ',num2str(mean(RMSE_elitist(:)))])
colorbar
subplot(2,4,7)
imshow(RMSE_fractional_im,[],'colormap', jet)
title(['RMSE fractional, mean: ',num2str(mean(RMSE_fractional(:)))])
colorbar
% subplot(2,4,8)
% imshow(RMSE_collaborative_im,[],'colormap', jet)
% title(['collaborative, mean: ',num2str(mean(RMSE_collaborative(:)))])
% colorbar
set(gcf,'color', 'white');

%%
% SAM (data)
SAM_FCLSU_im=reshape(SAM_FCLSU,m,n);
SAM_group_im=reshape(SAM_group,m,n);
SAM_elitist_im=reshape(SAM_elitist,m,n);
SAM_fractional_im=reshape(SAM_fractional,m,n);

figure
suptitle('The average SAM');
subplot(2,4,1)
imshow(SAM_g_im(:,:,2),[],'colormap',jet);
title(['SAM group(L0), mean: ',num2str(mean(SAM_g(:,2)))]);
colorbar
subplot(2,4,2)
imshow(SAM_e_im,[],'colormap',jet);
title(['SAM elitist(L0), mean: ',num2str(mean(SAM_e(:)))]);
colorbar
subplot(2,4,3)
imshow(SAM_mgs_im(:,:,1),[],'colormap',jet);
title(['SAM global(L0), mean: ',num2str(mean(SAM_mgs(:,1)))]);
colorbar
subplot(2,4,4)
imshow(SAM_FCLSU_im,[],'colormap', jet)
title(['SAM FCLSU, mean: ',num2str(mean(SAM_FCLSU(:)))])
colorbar
subplot(2,4,5)
imshow(SAM_group_im,[],'colormap', jet)
title(['SAM group, mean: ',num2str(mean(SAM_group(:)))])
colorbar
subplot(2,4,6)
imshow(SAM_elitist_im,[],'colormap', jet)
title(['SAM elitist, mean: ',num2str(mean(SAM_elitist(:)))])
colorbar
subplot(2,4,7)
imshow(SAM_fractional_im,[],'colormap', jet)
title(['SAM fractional, mean: ',num2str(mean(SAM_fractional(:)))])
colorbar
set(gcf,'color', 'white');

%%
% sort abundances for abundance matrix display as an image
minx = 7000;
maxx = 8000;

figure
% suptitle(['Abundance matrix of pixel ',num2str(minx),'to ',num2str(maxx)])
subplot(1,7,1)
imagesc(A_g(:,:,2))
xlim([minx,maxx])
title('Group(L0)','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)
subplot(1,7,2)
imagesc(A_e)
xlim([minx,maxx])
title('Elitist(L0)','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)
subplot(1,7,3)
imagesc(A_mgs(:,:,1))
xlim([minx,maxx])
title('Global(L0)','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)
subplot(1,7,4)
imagesc(A_FCLSU_bundle_sorted)
xlim([minx,maxx])
title('FCLSU bundles','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)
% subplot(1,6,6)
% imagesc(A_collaborative_sorted)
% xlim([minx,maxx])
% title('Collaborative','fontname','times','fontsize',12)
% set(gca,'fontname','times','fontsize',12)
subplot(1,7,5)
imagesc(A_group_sorted)
xlim([minx,maxx])
title('Group','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)
subplot(1,7,6)
imagesc(A_elitist_sorted)
xlim([minx,maxx])
title('Elitist','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)
subplot(1,7,7)
imagesc(A_fractional_sorted)
xlim([minx,maxx])
title('Fractional','fontname','times','fontsize',12)
set(gca,'fontname','times','fontsize',12)
set(gcf,'color','white')

%%
% Abundance maps according to groups
figure
suptitle('Abundance maps')
for p = 1:P
    
    subaxis(7,P,P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_g_im(:,:,p,2),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Group(L0)','fontname','times','fontsize',12)
    end
    colormap jet
    
    subaxis(7,P,2*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_e_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Elitist(L0)','fontname','times','fontsize',12)
    end
    colormap jet
    
    subaxis(7,P,3*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_mgs_im(:,:,p,1),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Global(L0)','fontname','times','fontsize',12)
    end
    colormap jet
    
    subaxis(7,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_FCLSU_bundle_corrected_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('FCLSU','fontname','times','fontsize',12)
    end
    colormap jet
    
%     subaxis(5,P,P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
%     imshow(A_collaborative_im(:,:,p),[],'colormap', jet)
%     set(gca,'clim',[0,1])
%     if p == 1
%         ylabel('Collaborative','fontname','times','fontsize',12)
%     end
%     colormap jet
    
    subaxis(7,P,4*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_group_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Group','fontname','times','fontsize',12)
    end
    colormap jet
    
    subaxis(7,P,5*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_elitist_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Elitist','fontname','times','fontsize',12)
    end
    colormap jet
    
    subaxis(7,P,6*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
    imshow(A_fractional_im(:,:,p),[],'colormap', jet)
    set(gca,'clim',[0,1])
    if p == 1
        ylabel('Fractional','fontname','times','fontsize',12)
        xlabel('Concrete','fontname','times','fontsize',12)
    elseif p == 2
        xlabel('Red Roofs','fontname','times','fontsize',12)
    elseif p == 3
        xlabel('Asphalt','fontname','times','fontsize',12)
    elseif p == 4
        xlabel('Vegetation','fontname','times','fontsize',12)
    else
        xlabel('Colored Structures','fontname','times','fontsize',12)
    end 
    
end
set(gcf,'color', 'white')


%%
% Abundance maps according to endmembers
for  i = 1:P
    
    EMs_cons = find(groups == i);
%     EM_cons = 1:bundle_nbr:P*bundle_nbr;
    EMs_L0_cons = find(spec_group == i);
    
    figure,
    suptitle(['Abundance maps of group ',num2str(i)])
    for p = 1:length(EMs_cons)
        
        subaxis(7,length(EMs_cons),p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
        imshow(A_FCLSU_bundle_corrected_full_im(:,:,EMs_cons(p)),[],'colormap',jet)
        set(gca,'clim',[0,1])
        if p == 1
            ylabel('FCLSU','fontname','times','fontsize',12)
        end
        
        subaxis(7,length(EMs_cons),length(EMs_cons)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
        imshow(A_g_full_im(:,:,EMs_L0_cons(p),2),[],'colormap',jet)
        set(gca,'clim',[0,1])
        if p == 1
            ylabel('Group(L0)','fontname','times','fontsize',12)
        end
        
        subaxis(7,length(EMs_cons),2*length(EMs_cons)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
        imshow(A_e_full_im(:,:,EMs_L0_cons(p)),[],'colormap',jet)
        set(gca,'clim',[0,1])
        if p == 1
            ylabel('Elitist(L0)','fontname','times','fontsize',12)
        end
        
        subaxis(7,length(EMs_cons),3*length(EMs_cons)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
        imshow(A_mgs_full_im(:,:,EMs_L0_cons(p),1),[],'colormap',jet)
        set(gca,'clim',[0,1])
        if p == 1
            ylabel('Global(L0)','fontname','times','fontsize',12)
        end
        
%         subaxis(5,length(EMs_cons),length(EMs_cons)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
%         imshow(A_collaborative_full_im(:,:,EMs_cons(p)),[],'colormap',jet)
%         set(gca,'clim',[0,1])
%         if p == 1
%             ylabel('Collaborative','fontname','times','fontsize',10)
%         end
        
        subaxis(7,length(EMs_cons),4*length(EMs_cons)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
        imshow(A_group_full_im(:,:,EMs_cons(p)),[],'colormap',jet)
        set(gca,'clim',[0,1])
        if p == 1
            ylabel('Group','fontname','times','fontsize',12)
        end
        
        
        subaxis(7,length(EMs_cons),5*length(EMs_cons)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
        imshow(A_elitist_full_im(:,:,EMs_cons(p)),[],'colormap',jet)
        set(gca,'clim',[0,1])
        if p == 1
            ylabel('Elitist','fontname','times','fontsize',12)
        end
        
        
        subaxis(7,length(EMs_cons),6*length(EMs_cons)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
        imshow(A_fractional_full_im(:,:,EMs_cons(p)),[],'colormap',jet)
        set(gca,'clim',[0,1])
        if p == 1
            ylabel('Fractional','fontname','times','fontsize',12)
        end
        set(gcf,'colormap', jet)
    end
    set(gcf,'color', 'white')
end

%% Reshape estimate
% rowDist=m*ones(1,n);
% est_sep=mat2cell(est_g',rowDist);
% est_data=cat(3,est_sep{:});
% est_data=permute(est_data,[1,3,2]);
% 
% figure, imshow(rescale(est_data(:,:,[57,30,20]),1),[])



