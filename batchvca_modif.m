function [ groups, components] = batchvca_modif( pixels, p, bundles, percent, clustering )

% pixels: data pixels, L*N
% p: number of endmember classes
% bundles: number of subsets of the data
% percent: percentage of the data taken each time
% clustering: either 'kmeans' or 'spectral_clustering'
% sampling is with replacement for now

    m = percent/100 * size(pixels,2);
    runs = 1;
    pixels_update = pixels;
    
%     B = [];  % change this for sampling without replacement
%     for b = 1:bundles
%         B = [B, vca_orig(pixels(:,randperm(size(pixels, 2), floor(m))), runs, p,false)];
%     end

    B = [];  % change this for sampling without replacement
    for b = 1:bundles
        [C,I] = datasample(pixels_update,floor(m),2,'Replace',false);
        B = [B, vca_orig(C, runs, p,false)];
        pixels_update(:,I) = [];
    end

    components = [];
    for i = 1:bundles
        components = [components, B(i).E];
    end
    
    % clustering part
    
    if strcmp(clustering,'kmeans')
        
        groups = kmeans(components',p,'distance','cosine');
        
    elseif strcmp(clustering,'spectral_clustering')
        
        Neighbors = 10;
        sigma = 1;
%         type = 1; % unnormalized
 type = 2; % normalized
        
        % to fill out
        % now for the clustering
        fprintf('Creating Similarity Graph...\n');
        SimGraph = SimGraph_NearestNeighbors(components, Neighbors, 1,sigma);
        
        try
            comps = graphconncomp(SimGraph, 'Directed', false);
            fprintf('- %d connected components found\n', comps);
        end
        
        fprintf('Clustering Data...\n');
        [~,groups]= SpectralClustering(SimGraph, p, type);
       
        elseif strcmp(clustering,'nystrom')
       
            alpha = 1;
            groups = nystromClustering(components',alpha,p);
            
    else
        error('Clustering method not recognized\n')
    end
    
%     comp_normed = components./repmat(sqrt(sum(components.^2,1)),size(components,1),1);
%     [U,~,~] = svd(comp_normed'*comp_normed);
% 
%     [groups, ~] = kmeans(U(:,[2:(p+1)]), p);
end
