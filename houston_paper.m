close all
clear 
dbstop if error

seed = 1; % super duper important

rng(seed) % super duper important


addpath(genpath('/home/administrateur/Documents/thèse/these_mai_2016/code/export_fig/export_fig-master/export_fig-master/'))
addpath(genpath('/home/administrateur/Documents/thèse/these_mai_2016/code/test EM var/joint unmixing of spatial datasets/synthesis/'))
addpath(genpath('/home/administrateur/Documents/thèse/these_mai_2016/code/test EM var/joint unmixing of spatial datasets/genabund/'))
addpath(genpath('/home/administrateur/Documents/thèse/these_mai_2016/code/data/pavia/'))
addpath(genpath('/home/administrateur/Documents/thèse/these_mai_2016/code/codes thèse/ID_algorithms/'))
addpath(genpath('/home/administrateur/Documents/thèse/these_mai_2016/code/codes thèse/synthesis/'))
addpath(genpath('/home/administrateur/Documents/thèse/these_mai_2016/code/japon/Manopt_4/manopt/'))
addpath(genpath('/home/administrateur/Documents/thèse/these_mai_2016/code/codes thèse/ELMM/ALS/toolbox/'))
addpath(genpath('/home/administrateur/Documents/these/these_mai_2016/code/export_fig/export_fig-master/export_fig-master/'))

load real_data_1

figure, imshow(rescale(data(:,:,[57,30,20]),1),[])

% data =rescale(x20100628_gipsa_area2(300:499,1:200,:),1) ;

EMs = figure;
imshow(rescale(data(:,:,[57,30,20]),1),[])
hold on

[m,n,L] = size(data);

%% Extract endmember candidates and cluster them into bundles

%  load bundle_test_paper_test.mat

X = reshape(data,m*n,L)';

imr= X;

P = 5;
% load bundle_test_paper.mat

bundle_nbr = 10;
percent = 10;
clustering = 'kmeans';
% clustering = 'spectral_clustering';
% clustering = 'nystrom'; % just to test function

% groups = kmeans(imr',P,'distance','cosine');
% bundle = imr;

[groups, bundle] = batchvca_modif(X, P, bundle_nbr, percent, clustering );
% [groups, bundle] = batchvca(X, P, bundle_nbr, percent);


pca_viz(imr,bundle)


% save('bundle_test_paper_test_no_sparse','groups','bundle')
% save('bundle_test_paper_test','groups','bundle')
% save('bundle_test_paper_test_difficult_synth_social','groups','bundle')
% compare visually using plots in each class + scatter plots
color_list{1} = [0 0 1];
color_list{2} = [0 1 0];
color_list{3} = [1 0 0];
color_list{4} = [1 0 1];
color_list{5} = [1 1 0];
color_list{6} = [0 1 1];
color_list{7} = [1 1 1];
color_list{8} = [0,0,0];
color_list{9} = [0.5,0.5,0];
color_list{10} = [0,0.5,0.5];
color_list{11} = [0.5,0,0.5];
color_list{12} = [0.75,0,0.25];
color_list{13} = [0.25,0,0.75];
color_list{14} = [0.5,0,0.75];
color_list{15} = [0.5,0.25,0.5];
% load bundle_test_paper_test_difficult_synth_social.mat
% load bundle_test_paper_test_no_sparse
% pca_viz_clusters(bundle,groups,color_list);
% pca_viz_clusters(reshape(S_true,L,N*P),ones(N*P,1),color_list);

% figure,
% for i = 1:P
% subplot(4,4,i)
% plot(squeeze(S_true(:,i,1)))
% end



% runs = 1;
% E = EIA_VCA(X,runs,P,false);
% vcaems = E.E;
% 
% figure 
% for i = 1:P
% subplot(4,5,i)
% plot(vcaems(:,i))
% ylim([0,1])
% end

%% unmix (empirical tuning of reg param)

% regular FCLSU on the bundles

tic
disp('FLCSU bundle')
A_FCLSU_bundle = FCLSU(X,bundle)';
t_FCLSU_bundle=toc;
A_init = A_FCLSU_bundle;

figure, imagesc(A_FCLSU_bundle)

% initialize params

%% RUN

% rho = 100;
rho = 10;
% tol_a = 10^(-10);
tol_a = 10^(-6);
fraction = 9/10;
% lambda = 100;

% maxiter_ADMM = 2000;
% maxiter_ADMM = 10000;
maxiter_ADMM = 1000;

tic
disp('FCLSU improved')
type = 'elitist';
[A_FCLSU_bundle_corrected, optim_struct_FCLSU] = social_unmixing(X,bundle,groups,A_init,0,rho,maxiter_ADMM,type,fraction,tol_a);
t_FCLSU_bundle_corrected=toc;

% lambda = 0.5;% OK

lambda = 2;

% maxiter_ADMM = 1000;
% group
tic
disp('group')
type = 'group';
rho=10;
[A_group, optim_struct_group] = social_unmixing(X,bundle,groups,A_init,lambda,rho,maxiter_ADMM,type,fraction,tol_a);
t_group=toc;

figure, plot(optim_struct_group.objective)
xlabel('group')

% elitist

% maxiter_ADMM = 1000;
% lambda =0.1;% NOT BAD
lambda =0.5;



tic
disp('elitist')
type = 'elitist';
[A_elitist, optim_struct_elitist] = social_unmixing(X,bundle,groups,A_init,lambda,rho,maxiter_ADMM,type,fraction,tol_a);
t_elitist=toc;
figure, plot(optim_struct_elitist.objective)
xlabel('elitist')
% fractional (use several values for p/q, start with 9/10.

% maxiter_ADMM = 1000;
fraction = 1/10;
% lambda = 0.1;
% lambda = 0.2; % NOT BAD

lambda = 0.4;

% lambda = 0.5; % SUPER SPARSE

tic
disp('fractional')
type = 'fractional';
[A_fractional, optim_struct_fractional] = social_unmixing(X,bundle,groups,A_init,lambda,rho,maxiter_ADMM,type,fraction,tol_a);
t_fractional=toc;


figure, plot(optim_struct_fractional.objective)
xlabel('fractional')


lambda = 1;
% tic
disp('collaborative')

%% Collaborative

type = 'asc';
tic
[A_collaborative] = ADMM_collaborative_unmixing(X,A_init,bundle,lambda,rho,maxiter_ADMM,type,tol_a,1);
t_collaborative=toc;
% sum the abundances within each class (bundle2global function)

[A_FCLSU_bundle_final,S_FCLSU_bundle] = bundle2global(A_FCLSU_bundle,bundle,groups);
[A_FCLSU_bundle_corrected_final,S_FCLSU_bundle_corrected] = bundle2global(A_FCLSU_bundle_corrected,bundle,groups);
[A_collaborative_final, S_collaborative] = bundle2global(A_collaborative,bundle,groups);
[A_group_final,S_group] = bundle2global(A_group,bundle,groups);
[A_elitist_final,S_elitist] = bundle2global(A_elitist,bundle,groups);
[A_fractional_final,S_fractional] = bundle2global(A_fractional,bundle,groups);


pca_viz_bundle(imr,reshape(S_FCLSU_bundle,L,m*n*P),bundle)

pca_viz_bundle(imr,reshape(S_group,L,m*n*P),bundle)

pca_viz_bundle(imr,reshape(S_elitist,L,m*n*P),bundle)

pca_viz_bundle(imr,reshape(S_fractional,L,m*n*P),bundle)

pca_viz_bundle(imr,reshape(S_collaborative,L,m*n*P),bundle)
% reconstruction errors

% H_FCLSU = S0*A_FCLSU; % reconstruction for FCLSU

H_FCLSU_bundle = bundle*A_FCLSU_bundle; % reconstruction for FCLSU
H_FCLSU_bundle_corrected = bundle*A_FCLSU_bundle_corrected; % reconstruction for FCLSU
H_group = bundle*A_group; % reconstruction for FCLSU
H_elitist = bundle*A_elitist; % reconstruction for FCLSU
H_fractional = bundle*A_fractional; % reconstruction for FCLSU
H_collaborative = bundle*A_collaborative;

%   RMSE_FCLSU = sqrt(1/L*sum((H_FCLSU'-data_r).^2,2));       
%   RMSE_FCLSU_im = reshape(RMSE_FCLSU,m,n);
  RMSE_FCLSU_bundle = sqrt(1/L*sum((H_FCLSU_bundle-X).^2,1));     
  RMSE_FCLSU_bundle_corrected = sqrt(1/L*sum((H_FCLSU_bundle_corrected-X).^2,1));       
  RMSE_FCLSU_bundle_im = reshape(RMSE_FCLSU_bundle,m,n);
  RMSE_FCLSU_bundle_corrected_im = reshape(RMSE_FCLSU_bundle,m,n);
  RMSE_group = sqrt(1/L*sum((H_group-X).^2,1));       
  RMSE_group_im = reshape(RMSE_group,m,n);
  RMSE_elitist = sqrt(1/L*sum((H_elitist-X).^2,1));       
  RMSE_elitist_im = reshape(RMSE_elitist,m,n);
  RMSE_fractional = sqrt(1/L*sum((H_fractional-X).^2,1));       
  RMSE_fractional_im = reshape(RMSE_fractional,m,n);
  RMSE_collaborative =  sqrt(1/L*sum((H_collaborative-X).^2,1));     
  RMSE_collaborative_im = reshape(RMSE_collaborative,m,n);
  
  N = m*n;
  
  SAM_FCLSU = zeros(N,1);
  SAM_group = zeros(N,1);
  SAM_elitist = zeros(N,1);
  SAM_fractional = zeros(N,1);
  SAM_collaborative = zeros(N,1);
  
for k = 1:N
        SAM_FCLSU(k) = 180/pi*real(acos((X(:,k)'*H_FCLSU_bundle_corrected(:,k))...
/(norm(X(:,k))*norm(H_FCLSU_bundle_corrected(:,k)))));
end

for k = 1:N
        SAM_group(k) = 180/pi*real(acos((X(:,k)'*H_group(:,k))...
/(norm(X(:,k))*norm(H_group(:,k)))));
end

for k = 1:N
        SAM_elitist(k) = 180/pi*real(acos((X(:,k)'*H_elitist(:,k))...
/(norm(X(:,k))*norm(H_elitist(:,k)))));
end

for k = 1:N
        SAM_fractional(k) = 180/pi*real(acos((X(:,k)'*H_fractional(:,k))...
/(norm(X(:,k))*norm(H_fractional(:,k)))));
end

for k = 1:N
        SAM_collaborative(k) = 180/pi*real(acos((X(:,k)'*H_collaborative(:,k))...
/(norm(X(:,k))*norm(H_fractional(:,k)))));
end


mean(SAM_FCLSU(:))
 mean(SAM_collaborative(:)) 
mean(SAM_group(:))
mean(SAM_elitist(:))
mean(SAM_fractional(:))
    

  
%% Qualitative results

% show RMSE (data)

    figure, 
    subplot(2,3,1)
  imshow(RMSE_FCLSU_bundle_im,[],'colormap', jet)
  title(['RMSE FCLSU bundle, mean: ',num2str(mean(RMSE_FCLSU_bundle(:)))])
  colorbar
      subplot(2,3,2)
  imshow(RMSE_collaborative_im,[],'colormap', jet)
  title(['collaborative, mean: ',num2str(mean(RMSE_collaborative(:)))])
  colorbar
      subplot(2,3,3)
  imshow(RMSE_group_im,[],'colormap', jet)
  title(['RMSE group, mean: ',num2str(mean(RMSE_group(:)))])
  colorbar
      subplot(2,3,4)
  imshow(RMSE_elitist_im,[],'colormap', jet)
  title(['RMSE elitist, mean: ',num2str(mean(RMSE_elitist(:)))])
  colorbar
      subplot(2,3,5)
  imshow(RMSE_fractional_im,[],'colormap', jet)
  title(['RMSE fractional, mean: ',num2str(mean(RMSE_fractional(:)))])
  colorbar
  set(gcf,'color', 'white')

% sort abundances for abundance matrix display as an image
A_FLCSU_bundle_sorted = zeros(size(A_FCLSU_bundle));
A_FLCSU_bundle_corrected_sorted = zeros(size(A_FCLSU_bundle));
A_group_sorted = zeros(size(A_group));
A_elitist_sorted = zeros(size(A_elitist));
A_fractional_sorted = zeros(size(A_fractional));
A_collaborative_sorted = zeros(size(A_collaborative));

group_sums = zeros(max(groups),1);

for i = 1:max(groups)
group_sums(i) = sum(groups == i);
end

group_cumsums = cumsum(group_sums);

for i = 1:max(groups)
    if i == 1
   A_FCLSU_bundle_sorted(1:group_cumsums(1),:) = A_FCLSU_bundle(groups == i,:); 
   A_FCLSU_bundle_corrected_sorted(1:group_cumsums(1),:) = A_FCLSU_bundle_corrected(groups == i,:); 
   A_group_sorted(1:group_cumsums(1),:) = A_group(groups == i,:); 
   A_elitist_sorted(1:group_cumsums(1),:) = A_elitist(groups == i,:); 
   A_fractional_sorted(1:group_cumsums(1),:) = A_fractional(groups == i,:); 
    A_collaborative_sorted(1:group_cumsums(1),:) = A_collaborative(groups == i,:); 
    else
   A_FCLSU_bundle_sorted(group_cumsums(i-1)+1:group_cumsums(i),:) = A_FCLSU_bundle(groups == i,:); 
   A_FCLSU_bundle_corrected_sorted(group_cumsums(i-1)+1:group_cumsums(i),:) = A_FCLSU_bundle_corrected(groups == i,:); 
   A_group_sorted(group_cumsums(i-1)+1:group_cumsums(i),:) = A_group(groups == i,:);  
   A_elitist_sorted(group_cumsums(i-1)+1:group_cumsums(i),:) = A_elitist(groups == i,:);  
   A_fractional_sorted(group_cumsums(i-1)+1:group_cumsums(i),:) = A_fractional(groups == i,:);  
    A_collaborative_sorted(group_cumsums(i-1)+1:group_cumsums(i),:) = A_collaborative(groups == i,:);  
    end
end

minx = 7000;
maxx = 8000;

figure,
subplot(1,5,1)
imagesc(A_FCLSU_bundle_sorted)
xlim([minx,maxx])
title('FCLSU bundles','fontname','times','fontsize',20)
set(gca,'fontname','times','fontsize',20)
subplot(1,5,2)
imagesc(A_collaborative_sorted)
xlim([minx,maxx])
title('Collaborative','fontname','times','fontsize',20)
set(gca,'fontname','times','fontsize',20)
subplot(1,5,3)
imagesc(A_group_sorted)
xlim([minx,maxx])
title('Group','fontname','times','fontsize',20)
set(gca,'fontname','times','fontsize',20)
subplot(1,5,4)
imagesc(A_elitist_sorted)
xlim([minx,maxx])
title('Elitist','fontname','times','fontsize',20)
set(gca,'fontname','times','fontsize',20)
subplot(1,5,5)
imagesc(A_fractional_sorted)
xlim([minx,maxx])
title('Fractional','fontname','times','fontsize',20)
set(gca,'fontname','times','fontsize',20)
set(gcf,'color','white')




A_FCLSU_bundle_im = reshape(A_FCLSU_bundle_final',m,n,P);
A_FCLSU_bundle_corrected_im = reshape(A_FCLSU_bundle_corrected_final',m,n,P);
A_group_im = reshape(A_group_final',m,n,P);
A_elitist_im = reshape(A_elitist_final',m,n,P);
A_fractional_im = reshape(A_fractional_final',m,n,P);
A_collaborative_im =reshape(A_collaborative_final',m,n,P);

figure,
for p = 1:P
   
       subaxis(5,P,p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_FCLSU_bundle_corrected_im(:,:,p),[],'colormap', jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('FCLSU','fontname','times','fontsize',15)
   end
   colormap jet

       subaxis(5,P,P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_collaborative_im(:,:,p),[],'colormap', jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('Collaborative','fontname','times','fontsize',15)
   end
   colormap jet
   
    subaxis(5,P,2*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_group_im(:,:,p),[],'colormap', jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('Group','fontname','times','fontsize',15)
   end
   colormap jet

     subaxis(5,P,3*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_elitist_im(:,:,p),[],'colormap', jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('Elitist','fontname','times','fontsize',15)
   end
   colormap jet
   
    subaxis(5,P,4*P+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_fractional_im(:,:,p),[],'colormap', jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('Fractional','fontname','times','fontsize',15)
   xlabel('Concrete','fontname','times','fontsize',15)
   elseif p == 2
       xlabel('Red Roofs','fontname','times','fontsize',15)
   elseif p == 3
       xlabel('Vegetation','fontname','times','fontsize',15)
   elseif p == 4
       xlabel('Asphalt','fontname','times','fontsize',15)
   else 
       xlabel('Colored Structures','fontname','times','fontsize',15)
   end
   
   
  
end
 set(gcf,'color', 'white')

%%

A_FCLSU_bundle_full_im = reshape(A_FCLSU_bundle',m,n,P*bundle_nbr);
A_FCLSU_bundle_corrected_full_im = reshape(A_FCLSU_bundle_corrected',m,n,P*bundle_nbr);
A_group_full_im = reshape(A_group',m,n,P*bundle_nbr);
A_elitist_full_im = reshape(A_elitist',m,n,P*bundle_nbr);
A_fractional_full_im = reshape(A_fractional',m,n,P*bundle_nbr);
A_collaborative_full_im = reshape(A_collaborative',m,n,P*bundle_nbr);

% EMs_cons = [4,15,25,35,40,48];

for  i = 1:P

EMs_cons = find(groups == i);
% EM_cons = 1:bundle_nbr:P*bundle_nbr;


figure,
for p = 1:length(EMs_cons)

       subaxis(5,length(EMs_cons),p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_FCLSU_bundle_corrected_full_im(:,:,EMs_cons(p)),[],'colormap',jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('FCLSU','fontname','times','fontsize',10)
   end
   

       subaxis(5,length(EMs_cons),length(EMs_cons)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_collaborative_full_im(:,:,EMs_cons(p)),[],'colormap',jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('Collaborative','fontname','times','fontsize',10)
   end
   
    subaxis(5,length(EMs_cons),2*length(EMs_cons)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_group_full_im(:,:,EMs_cons(p)),[],'colormap',jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('Group','fontname','times','fontsize',10)
   end
   

     subaxis(5,length(EMs_cons),3*length(EMs_cons)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_elitist_full_im(:,:,EMs_cons(p)),[],'colormap',jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('Elitist','fontname','times','fontsize',10)
   end
   
   
    subaxis(5,length(EMs_cons),4*length(EMs_cons)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_fractional_full_im(:,:,EMs_cons(p)),[],'colormap',jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('Fractional','fontname','times','fontsize',10)
   end
  set(gcf,'colormap', jet)
end
 set(gcf,'color', 'white')
end

%%


EMs_cons = find(groups == 4);
% EM_cons = 1:bundle_nbr:P*bundle_nbr;
ems_sel = [1,3,6,7,9];

figure,
for p = 1:length(ems_sel)


       subaxis(5,length(ems_sel),p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_FCLSU_bundle_corrected_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('FCLSU','fontname','times','fontsize',15)
   end
   
       subaxis(5,length(ems_sel),length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_collaborative_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('Collaborative','fontname','times','fontsize',15)
   end
   
    subaxis(5,length(ems_sel),2*length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_group_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('Group','fontname','times','fontsize',15)
   end
   

     subaxis(5,length(ems_sel),3*length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_elitist_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('Elitist','fontname','times','fontsize',15)
   end
   
   
    subaxis(5,length(ems_sel),4*length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_fractional_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('Fractional','fontname','times','fontsize',15)
   end
  set(gcf,'colormap', jet)
end
 set(gcf,'color', 'white')
 
 
 %%
 EMs_cons = find(groups == 2);
% EM_cons = 1:bundle_nbr:P*bundle_nbr;
ems_sel = [1,3,5,9,4];

figure,
for p = 1:length(ems_sel)


       subaxis(5,length(ems_sel),p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_FCLSU_bundle_corrected_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('FCLSU','fontname','times','fontsize',15)
   end
   
       subaxis(5,length(ems_sel),length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_collaborative_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('Collaborative','fontname','times','fontsize',15)
   end
   
    subaxis(5,length(ems_sel),2*length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_group_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('Group','fontname','times','fontsize',15)
   end
   

     subaxis(5,length(ems_sel),3*length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_elitist_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('Elitist','fontname','times','fontsize',15)
   end
   
   
    subaxis(5,length(ems_sel),4*length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_fractional_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('Fractional','fontname','times','fontsize',15)
   end
  set(gcf,'colormap', jet)
end
 set(gcf,'color', 'white')
 
 %%
 
 EMs_cons = find(groups == 1);
% EM_cons = 1:bundle_nbr:P*bundle_nbr;
ems_sel = [1,4,6,10,11];

figure,
for p = 1:length(ems_sel)


       subaxis(5,length(ems_sel),p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_FCLSU_bundle_corrected_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('FCLSU','fontname','times','fontsize',15)
   end
   
       subaxis(5,length(ems_sel),length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_collaborative_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('Collaborative','fontname','times','fontsize',15)
   end
   
    subaxis(5,length(ems_sel),2*length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_group_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('Group','fontname','times','fontsize',15)
   end
   

     subaxis(5,length(ems_sel),3*length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_elitist_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('Elitist','fontname','times','fontsize',15)
   end
   
   
    subaxis(5,length(ems_sel),4*length(ems_sel)+p,'SpacingVertical',0.01,'SpacingHorizontal',0.01)
   imshow(A_fractional_full_im(:,:,EMs_cons(ems_sel(p))),[],'colormap',jet)
   set(gca,'clim',[0,1])
   if p == 1
   ylabel('Fractional','fontname','times','fontsize',15)
   end
  set(gcf,'colormap', jet)
end
 set(gcf,'color', 'white')

 %%
 
%  idx_display = randi(N,1000);
%  
%  figure 
% for i = 1:P
% subplot(2,5,i)
% plot(bundle(:,groups == i))
% ylim([0,1])
% xlim([1,144])
% subplot(2,5,i)
% plot(squeeze(S_fractional(:,i,idx_display)))
% end

%%

pca_viz_global(imr,reshape(S_FCLSU_bundle,L,P*N),reshape(S_collaborative,L,P*N),reshape(S_group,L,P*N),reshape(S_elitist,L,P*N),reshape(S_fractional,L,P*N),bundle)

 
 

