% REEFMOD-GBR model run script
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, 02/2021
%__________________________________________________________________________

function [META, REEF, RESULT, RECORD] = f_multiple_reef(OPTIONS, RESTORATION, nb_time_steps, simul)

PARAMETERS

load('GBR_REEF_POLYGONS.mat')
META.reef_ID = [1:3806]; % ID of the modelled reefs
META.area_habitat = GBR_REEFS.HabitatAreaKm2(META.reef_ID);
META.nb_reefs = length(META.reef_ID);
META.nb_time_steps = nb_time_steps;

REEF.herbivory = 1; % full grazing
ALGAL.nb_step_algal_dynamics = 1 ; %%%%% ONLY TO SPEED-UP THE CODE WITH FULL GRAZING  (otherwise use 6)

%% Disturbances
META.doing_bleaching = OPTIONS.doing_bleaching ;
META.deterministic_bleaching = 0; % option for generating deterministic (1) or random (0) mortalites from DHWs
META.DHW_threshold = 3 ; % DHW threshold that triggers bleaching

% Elevate bleaching to align with the shallow (-2m) mortality as recorded by Hughes et al. (2018)
CORAL.sensitivity_bleaching = 2*CORAL.sensitivity_bleaching; % shallow bleaching
% CORAL.sensitivity_bleaching = CORAL.sensitivity_bleaching; % deep bleaching

META.doing_hurricanes = OPTIONS.doing_cyclones ;
META.random_hurricanes = 0; % put 0 to apply a specific scenario (1 imposes random occurence)
META.deterministic_hurricane_mortality = 0; % option for generating deterministic (1) or random (0) mortalites from cyclone cat

META.doing_COTS = OPTIONS.doing_COTS ;
META.doing_COTS_control = 0 ;
REEF_COTS =[]; % (required if not doing CoTS)
META.randomize_initial_COTS_densities = 2 ; % Distri for random generation of COTS densities (based on CoCoNet averages)
% 1 for Gaussian, 2 for Poisson

META.doing_water_quality = OPTIONS.doing_WQ;
META.randomize_WQ_chronology = 0;

META.doing_size_frequency = OPTIONS.doing_size_frequency; % Track coral colonies sizes (SLOW!)

%% Connectivity
META.doing_coral_connectivity = 1 ;

META.recruitment_type = 1; % turn into 0 for fixed recruitment (but then connect and genetics won't work)
META.coral_immigration = 1e6*ones(1,META.nb_time_steps) ; % Forced larval input for a 400m2 area - works only if connectivity is OFF

% Parameter a of the B-H function (same for all reefs)
CORAL.BH_alpha = 15*CORAL.prop_settlers;
CORAL.BH_beta = 5*1e6*ones(1,6); % calibrated with with META.coral_min_selfseed = 0.28

% Force self-seeding of coral larvae
META.coral_min_selfseed = 0.28 ; % relative proportion of produced larvae that stay on the reef (Helix experiment)

%% Rubble
META.tracking_rubble = 1;
META.rubble_decay_rate = 0.128 ; % 2/3 stabilized after 4 years
META.convert_rubble = 1; % Conversion factor from coral loss to rubble cover

%% Restoration options
META.nb_restored_reefs = RESTORATION.nb_reefs_restored ;

if META.nb_restored_reefs>0   
    settings_RESTORATION
end

if RESTORATION.doing_cooling==0
    RESTORATION.cooling=0; % otherwise requires cooling factor for each time step
end

%% Genetic adaptation
META.doing_genetics = OPTIONS.doing_adaptation ;

if META.doing_genetics==1
    
    settings_GENETICS;
    
    META.genetics.SIGMA_COLD = OPTIONS.adaptation_parms(1) + [ 0 0 0 0 0 0 ]; % for the cold side (when temp<Topt)
    META.genetics.SIGMA_HOT = OPTIONS.adaptation_parms(2) + [ 0 0 0 0 0 0 ]; % for the hot side (when temp>Topt)
    META.genetics.esd = OPTIONS.adaptation_parms(3) + [ 0 0 0 0 0 0 ]; % SD of environmental effect on fitness (mean=0), on top of genetics.
    META.genetics.enhanced_tolerance = OPTIONS.thermal_tolerance ;
    
    % Load the pre-adapted pool of QTL
    load(['QTL_pool_sigma_c' num2str(OPTIONS.adaptation_parms(1)) '_h' num2str(OPTIONS.adaptation_parms(2))...
        '_esd' num2str(OPTIONS.adaptation_parms(3)) '.mat']);
    
    CORAL.growth_rate = CORAL.growth_rate/mean_fitness; % average fitness across the GBR after burnin
    % this allows adjusting fitness so that mean individual growth rate is close to empirical values from literature
    % Note this average fitness is relative to the local environment corals are adapted to, ie we are not 
    % modelling latitudinal differences in growth rates (100% fitness in the far South gives the same growth rate than 
    % in the far North, while in reality coral growth declines at with latitude
    
end

%% INITIALISATION
% rng('shuffle')
rng(simul) % to get a repeatable scheme of random number generation in RAND, RANDI, RANDN

INITIALISATION

settings_GBR

MULTIPLE_REEF_SETUP


clear QTL1000 CYCLONE_CAT DHW GBR_REEF_POP SST_GBR SST_baseline GBR_REEFS MMM_CoRTAD
clear init_coral_cover init_rubble_cover init_sand_cover reef n id1 id2

[RESULT, RECORD] = f_runmodel(META, REEF, CORAL, ALGAL, CONNECT, REEF_POP, REEF_COTS) ;
