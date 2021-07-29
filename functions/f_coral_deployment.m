% -------------------------------------------------------------------------
% Y.-M. Bozec, MSEL, created Aug 2018.
% Deployment of corals for restoration (with genetics)
% -------------------------------------------------------------------------

function   [coral, algal, genes, ID_colony_tracking, total_transplanted] = ...
    f_coral_deployment(coral, algal, genes, REEF, ID_colony_tracking, META)

select_species = find(META.restore_coral_per_cell ~= 0);

list_cell = find(REEF.grazable_cell==1);
r = randperm(length(list_cell)) ;
list_cell2 = list_cell(r) ;

avail_space_cm2 = algal(1).cover_cm2 + algal(4).cover_cm2 ; % can be only outplanted on EAM and TURF

total_transplanted = zeros(size(coral));

for s = select_species
    
    [l,c] = size(coral(s).cover_cm2);   
    id1 = zeros(l,c);
    id1(coral(s).cover_cm2>0)=1;
    colony_count = sum(id1,2) ; % number of colonies of a given species for each cell
    
    transplant_layer = zeros(l, META.max_colonies-c);
    transplant_size_cm2 = floor(pi*(META.restore_coral_diameter(s)/2)^2);
    nb_transplants = zeros(l,1);      
            
    if META.restore_random_density == 1
        
        nb_transplants(list_cell2,1) = poissrnd(META.restore_coral_per_cell(s),1,length(list_cell2)) ;
        nb_transplants(nb_transplants > 5) = 5;
        
        check = META.max_colonies - nb_transplants - colony_count ;
        nb_transplants(check<0) = META.max_colonies ;

    else
        
        nb_transplants(list_cell2,1) = floor(META.restore_coral_per_cell(s))*ones(length(list_cell2),1) ;
        % should first determine total number of transplants over the grid then distribute
        
    end
      
    for j = 1:size(transplant_layer,2)
        
        % Need to implement cap of coral cover regarding available space
        do_transplant = zeros(l,1);
        do_transplant(nb_transplants>0)=1;
        do_transplant(avail_space_cm2(list_cell2) < transplant_size_cm2) = 0 ; % cannot exceed available space
        transplant_layer(list_cell2,j) = do_transplant(list_cell2)*transplant_size_cm2 ;
        nb_transplants = nb_transplants - do_transplant ;
        
    end

    cols = sum(transplant_layer,1);
    transplant_layer = transplant_layer(:,cols~=0) ;
    total_transplanted(s) = sum(cols/transplant_size_cm2);

    transplant_ID = zeros(size(transplant_layer));
    transplant_ID(transplant_layer>0)=[1:1:total_transplanted(1,s)]+ID_colony_tracking(s);

    coral(s).cover_cm2 = [coral(s).cover_cm2 transplant_layer(:,cols~=0)];
    coral(s).colony_ID = [coral(s).colony_ID transplant_ID]; % marked with -1
    
    if META.doing_clades == 1        
        coral(s).clade = [coral(s).clade transplant_layer/transplant_size_cm2]; %transplant clade 1        
    end
    
    ID_colony_tracking(1,s) = ID_colony_tracking(1,s)+total_transplanted(1,s);% Record the new ID max

%     algal(4).cover_cm2 = algal(4).cover_cm2 - sum(transplant_layer,2); % should overgrow EAM as well
    algal(1).cover_cm2 = algal(1).cover_cm2 - sum(transplant_layer,2); % should overgrow EAM as well
    
    % Assign TQLs
    if META.doing_genetics == 1
         
%         enhance_tolerance = META.genetics.enhanced_tolerance(s)/(size(META.genetics.QTL_pool,2)*size(META.genetics.QTL_pool,3));
        enhance_tolerance = META.genetics.enhanced_tolerance(s)/(size(REEF.coral(s).QTL_pool_IN,2)*size(REEF.coral(s).QTL_pool_IN,3));
        
%         idx = randi(size(META.genetics.QTL_pool,1),total_transplanted(1,s),1);
%         TEMP_QTL = META.genetics.QTL_pool(idx,:,:) + enhance_tolerance ;
%         
%         genes(s).QTLs = [genes(s).QTLs ; TEMP_QTL] ;
%         genes(s).list_coral_ID = [genes(s).list_coral_ID ; transplant_ID(transplant_ID~=0)];
%         
%         % Compute genotypes with heritability
%         % First calculate breeding value
%         BREED = sum(sum(TEMP_QTL,3),2) ;
%         % Then deeermine phenotype following heritability (esd=0 implies perfect heritability)
%         tmp_phenotypes = BREED + normrnd(0, META.genetics.esd(s), size(BREED)) ;
%         genes(s).phenotypes = [genes(s).phenotypes ; tmp_phenotypes] ;
%         
        
        %%%%%% NEW WAY
        idx = randi(size(REEF.coral(s).QTL_pool_IN,1),total_transplanted(1,s),1);
        TEMP_QTL = REEF.coral(s).QTL_pool_IN(idx,:,:) + enhance_tolerance ;
        
        genes(s).QTLs = [genes(s).QTLs ; TEMP_QTL] ;
        genes(s).list_coral_ID = [genes(s).list_coral_ID ; transplant_ID(transplant_ID~=0)];
        
        % Compute genotypes with heritability
        % First calculate breeding value
        BREED = sum(sum(TEMP_QTL,3),2) ;
        % Then determine phenotype following heritability (esd=0 implies perfect heritability)
        Env_effect = normrnd(0, META.genetics.esd(s), size(BREED));
%         Env_effect(Env_effect>4) = 4; % Limit the environmental effect (add no more that +6 deg C to BREED)
              
        % Then deeermine phenotype following heritability (esd=0 implies perfect heritability)
        tmp_phenotypes = BREED + Env_effect ;
        
        max_phenotype = floor(max(META.List_Topt)-REEF.Topt_baseline);
        tmp_phenotypes(tmp_phenotypes > max_phenotype) = max_phenotype ;        
       
        genes(s).phenotypes = [genes(s).phenotypes ; tmp_phenotypes] ;
        
    end
    
end
    