%__________________________________________________________________________
%
% REEFMOD-GBR MAIN SCRIPT
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, 02/2021
%__________________________________________________________________________
clear

NB_TIME_STEPS = 26; % HINDCAST: summer 2008 to winter 2020 (initialisation = winter 2007)
% NB_TIME_STEPS = 26+100; % FORECAST summer 2008 - winter 2070 (always re-run the hindcast before future projections

NB_SIMULATIONS = 40; % Number of repeated runs

OutputFilename = 'R0_HINDCAST_GBR.mat'; options = [ 1 1 1 1];
% OutputFilename = 'R1_HINDCAST_GBR_NO_CYCL.mat';   options = [ 0 1 1 1 ];
% OutputFilename = 'R2_HINDCAST_GBR_NO_BLEACH.mat'; options = [ 1 0 1 1 ];
% OutputFilename = 'R3_HINDCAST_GBR_NO_COTS.mat';   options = [ 1 1 0 1 ];
% OutputFilename = 'R4_HINDCAST_GBR_NO_WQ.mat';     options = [ 1 1 1 0 ];
SaveDir = 'outputs';

%% --------------------------------------------------------------------------------
% Options: yes(1)/no(0)
OPTIONS.doing_cyclones = options(1);
OPTIONS.doing_bleaching = options(2) ; 
OPTIONS.doing_COTS = options(3);
OPTIONS.doing_WQ = options(4);

OPTIONS.doing_size_frequency = 1; % for tracking population size structure (incl. juveniles)

OPTIONS.doing_adaptation = 0 ;
OPTIONS.adaptation_parms = [ 1 1 1.5 ]; % sigma cold/sigma hot/esd
% OPTIONS.adaptation_parms = [ 2 2 0.5 ]; % sigma cold/sigma hot/esd

RESTORATION.nb_reefs_restored = 0 ;
RESTORATION.doing_cooling = 0 ;
RESTORATION.thermal_tolerance = 0 ;


%% --------------------------------------------------------------------------------
OUTPUTS = struct('REEF', [],'RESULT', [],'RECORD', []);
TEMP_META = struct('META', []);

% parfor i=1:NB_SIMULATIONS
for i=1:NB_SIMULATIONS
    
    i
    
    [meta, REEF, RESULT, RECORD] = f_multiple_reef(OPTIONS, RESTORATION, NB_TIME_STEPS, i);
      
    OUTPUTS(i).RESULT = RESULT ;
    OUTPUTS(i).RECORD = RECORD ;
    OUTPUTS(i).REEF = REEF ;
    TEMP_META(i).META = meta ;
end

META = TEMP_META(1).META; % Keep only one META because common to all simulations
clear TEMP_META ADAPT i 
                
%% --------------------------------------------------------------------------------
%% Memory allocation for coral outputs
coral_cover_per_taxa = zeros(NB_SIMULATIONS,META.nb_reefs,META.nb_time_steps+1,META.nb_coral_types,'single');
% coral_outplanted = coral_cover_per_taxa;
coral_larval_supply = coral_cover_per_taxa;
coral_recruits = coral_cover_per_taxa;
coral_offsprings = coral_cover_per_taxa;
    
if OPTIONS.doing_size_frequency == 1   
    nb_juv = zeros(NB_SIMULATIONS, NB_REEFS, NB_TIME_STEPS+1, 6, 4) ;
    nb_adol = zeros(NB_SIMULATIONS,META.nb_reefs,META.nb_time_steps+1,META.nb_coral_types, 13,'single') ;
    nb_adult = zeros(NB_SIMULATIONS,META.nb_reefs,META.nb_time_steps+1,META.nb_coral_types, 9,'single') ; 
end

% Cover loss following stressors
coral_cover_lost_bleaching = zeros(NB_SIMULATIONS,META.nb_reefs,META.nb_time_steps,META.nb_coral_types,'single');
coral_cover_lost_cyclones = coral_cover_lost_bleaching;
coral_cover_lost_COTS = coral_cover_lost_bleaching;

% Other variables
rubble = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'single');
nongrazable = zeros(NB_SIMULATIONS, META.nb_reefs,'single');

%% Memory allocation for CoTS outputs
COTS_mantatow = zeros(NB_SIMULATIONS,META.nb_reefs,META.nb_time_steps+1,'single');
COTS_densities = zeros(NB_SIMULATIONS,META.nb_reefs,META.nb_time_steps+1,16,'single'); % 16 age classes
COTS_settler_density = COTS_mantatow;
COTS_larval_supply = COTS_mantatow;
COTS_larval_output = COTS_mantatow;

%% Populate outputs

for simul = 1:NB_SIMULATIONS
    
    coral_cover_per_taxa(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_pct2D));
%     coral_outplanted(simul,:,2:end,:) = squeeze(cat(4,OUTPUTS(simul).RECORD.total_transplanted));
    coral_larval_supply(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_larval_supply));
    coral_recruits(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_settler_count));
    coral_offsprings(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_total_fecundity));
    
    if do_size_frequency == 1
        nb_juv(simul,:,:,:,:)= squeeze(cat(4,OUTPUTS(simul).RESULT.coral_juv_count(:,:,:,:)));
        nb_adol(simul,:,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_adol_count(:,:,:,:)));
        nb_adult(simul,:,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_adult_count(:,:,:,:)));
    end

    % Assuming 0.6 CoTS per grid ~ 0.22 CoTS per tow 
    % (0.22 per tow is equivalent to 1500 COTS per km2, so that 1 COTS per grid (400m2) is equivalent to 0.22*2500/1500
    COTS_mantatow(simul,:,:) = (0.22/0.6)*squeeze(cat(4,OUTPUTS(simul).RESULT.COTS_total_perceived_density));
    COTS_densities(simul,:,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_all_densities); % Density for 400m2
    COTS_settler_density(simul,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_settler_densities); % Density for 400m2
    COTS_larval_supply(simul,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_larval_supply); % Density for 400m2
    COTS_larval_output(simul,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_larval_output); % Density for 400m2
    
    coral_cover_lost_bleaching(simul,:,:,:) = squeeze(OUTPUTS(simul).RECORD.coral_pct2D_lost_bleaching);
    coral_cover_lost_cyclones(simul,:,:,:) = squeeze(OUTPUTS(simul).RECORD.coral_pct2D_lost_cyclones);
    coral_cover_lost_COTS(simul,:,:,:) = squeeze(OUTPUTS(simul).RECORD.coral_pct2D_lost_COTS);

    rubble(simul,:,:)= squeeze(OUTPUTS(simul).RESULT.rubble_cover_pct2D); 
    nongrazable(simul,:) = squeeze(cat(4,OUTPUTS(simul).REEF.nongrazable_substratum));        
    
end

clear OUTPUTS simul s 

save([SaveDir OutputFilename])
