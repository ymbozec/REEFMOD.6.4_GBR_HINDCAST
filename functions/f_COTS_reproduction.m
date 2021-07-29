%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Sep 2017.
%
% This estimates the production of COTS larvae from a a given population size and age structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ egg_production ] = f_COTS_reproduction( RESULT, META, reef, time )


egg_production = sum(squeeze(RESULT.COTS_all_densities(reef,time,:)).*META.COTS_fecundity');

% Below is attempt at introducing more complex relationships between
% reproductive output and CoTS density following Babcock et al. 2014.
% Finally ignored here for the sake of simplicity, also because benefits of
% the code below is unclear: a power function with b<1 (here 0.61) leads to
% a plateau of fertilization rate with coTS density, but a plateau is also
% explicit in the formulation of the Beverton-Hold function linking CoTS
% recruits and larval supply (at sink). The difference is that reefs with
% high coTS densities are finally limited in the number of larvae produced.
% Other things to clarify: Allee effect mentionned in Babcock et al. 2014
% and exponential relationship between zygote production and CoTS density -
% clearly different dynamics. But these relationships are meant to reflect 
% the effect of CoTS aggregations which are not well captured by our model:
% 1 reef = 20x20m reef susbtratum, with CoTS density reflecting 6 month of 
% invasion.

% density_matures = sum(RESULT.COTS_all_densities(reef,time,5:end)) ;
% 
% fertilization_success = 0.14 * ((10^8)*density_matures/META.total_area_cm2)^0.61 ; % Babcock et al. (2014) with COTs density per hectare
% fertilization_success(fertilization_success>0.9)=0.9;  %fertilization cap (max obtained by Babcock)
% 
% egg_production = sum(squeeze(RESULT.COTS_all_densities(reef,time,:)).*META.COTS_fecundity')*fertilization_success;

