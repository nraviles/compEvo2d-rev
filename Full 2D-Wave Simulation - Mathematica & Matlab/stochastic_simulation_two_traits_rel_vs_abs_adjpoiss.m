function [pop, fit_a, fit_r] = stochastic_simulation_two_traits_rel_vs_abs_adjpoiss( ...
                                T,init_pop,init_fit_a,init_fit_r,b,d,sa,ua,uad, ...
                                sr,ur,urd,E,steps,collect_data,outputfile)
% The code below has been modified from the source code made
% availabe by Pearce MT and Fisher DS, obtained from:
% 
% https://datadryad.org/resource/doi:10.5061/dryad.f36v6 
%
% Stochastic simulation for two trait model, with relative vs absolute fitness. 
% Each step is one generation. Rates (such as s,u,r) are per generation. 
% The expected size of each subpop after selection, mutation, and mating is computed.
% If expected size < :cutoff: then the size is drawn from Poisson
% distribution with mean equal to the expected size. Otherwise the
% subpopulation size = expected size. 

% output :v, v1, v2:        total rate of adaptation, and rate of adaptation in traits 1 and 2. 
% output :varx, vary, cov:  time averaged variances and covariance in traits 1 and 2.

% input :T: available territories
% input :sa: effect size of beneficial mutation in trait 1 (absolute fitness), reduces the death rate of an individual. 
% input :ua: "increment" in mutation rate per locus of trait 1 (absolute fitness).
% input :uad: deleterious mutation rate per locus of trait 1 (absolute fitness).
% input :sr: effect size of beneficial mutation in trait 2 (relative fitness), increases the competitive ability of an individual. 
% input :ur: mutation rate per locus of trait 2 (relative fitness).
% input :urd: deleterious mutation rate per locus of trait 2 (relative fitness).
% input :E: rate of environmental change (1/E gives mean time between changes)
% input :steps: maximum number of steps for simulation.
% input :collect_data: true/false - collect detailed data on 2d distr. per generation
% input :start_time: start time for collecting detailed data on 2d distribution
% input :end_time: end time for collecting detailed data on 2d distribution
% input :outputfile: string with filename where detailed data will be stored
%
%                 rel fit 
%               ------------>
%    |    [ * * * * * * * * * * ]
%    |    [ * * * * * * * * * * ]
% abs fit [ * * * * * * * * * * ]
%    |    [ * * * * * * * * * * ]
%    |    [ * * * * * * * * * * ]
%    v    [ * * * * * * * * * * ]
%         [ * * * * * * * * * * ]
%
% NOTE: sa should be picked to ensure 1/(sa*d) is an integer, and in fact
% number of fitness classes  = 1/(sa*d).

%digits(16);

% initialize variables
pop = init_pop;                      % abundances of a classes
ext_time = 0;
lead_distr = [-round(1/(sa*d)):1:0];
lead_count = zeros(size([-round(1/(sa*d)):1:0]));

fit_a=init_fit_a;                    % trait 1 beneficial mutations (absolute fitness)
fit_r=init_fit_r;                    % trait 2 beneficial mutations (relative fitness)

% Absolute fitness is tracked by number of mutations from optimal genotype,
% i.e fit_a = [-4 .... -1 0]. Relative fitness is unbounded and will
% consist of arrays like fit_r = [10 11 .... 100].

cutoff = 1000/sr;%10/sr;       % population cutoff for stochasticity

if (collect_data)           % store parameters used in simulation
    fileID0 = fopen([outputfile '-0.txt'],'w');
    % fprintf(fileID0,'%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f',T,sa,ua,uad,sr,ur,urd,E,steps);
    fclose(fileID0);
    fileID1 = fopen([outputfile '-1.txt'],'w'); %file for extinction time versus mean time between env changes
    fileID2 = fopen([outputfile '-2.txt'],'w'); %file for classes
    fileID3 = fopen([outputfile '-3.txt'],'w'); %file for abundances
end

% Sample time of first ENVIRONMENTAL CHANGE
te = steps + 1; %max([round(exprnd(1/E)) 1]);

% Main loop for simulation of each generation
while ( (sum(sum(pop))>0) && (ext_time <= steps) )
    
    lead_count(lead_distr == max(fit_a)) = lead_count(lead_distr == max(fit_a))+1;
    
    % SELECTION STEP USING BERTRAM & MASEL LOTTERY MODEL
    % note that Wi, W_bar, wi are fitness for current set of classes, not
    % for the new set after mutations down below.
    
    [pop,delta_ni,Wi,W_bar,wi,c_bar] = BM_lottery_model_adjpoiss(pop,fit_a,fit_r,T,b,d,sa,sr); 
    
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    
    % MUTATIONS AND DRIFT STEP
    [pop,fit_a,fit_r] = two_dim_abs_vs_rel_mutations_adjpoiss(pop,delta_ni,fit_a,fit_r,ua,uad,ur,urd,cutoff);
    
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    
    % CHECK FOR NEXT ENVIRONMENTAL CHANGE, OR CONTINUE COUNTING DOWN
    if(te == 0)
        fit_a = fit_a - 1;          %absolute fitness
        te = max([round(exprnd(1/E)) 1]);
    else
        te = te -1;
    end
    
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------    
    
%     % burn in period calculations
%     if timestep > 3000
%         
%     else
%         if timestep==3000
% 
%         end
%     end
    
    if( collect_data && (sum(sum(pop))==0) )
        
        for i=1:size(lead_distr,2)
            % print data to output files, need: times,mean_fit,fit_var,fit_cov,pop_load,dcov_dt,vU_thry,v2U_thry
            % fprintf(fileID1,'%i,%i\n',lead_distr(i),lead_count(i));
        end
        
        for i=1:size(pop,1)
            for j=1:size(pop,2)
                if(pop(i,j)>0)
                    % fprintf(fileID2,'[%i,%i],',fit_a(i),fit_r(j));
                    % fprintf(fileID3,'%f,',pop(i,j));
                end
            end
        end
        
        % fprintf(fileID2,'\n');
        % fprintf(fileID3,'\n');
    end
    
    ext_time = ext_time+1;
    %ext_time
end

% close output files
if(collect_data)
    fclose(fileID1);
    fclose(fileID2);
    fclose(fileID3);
end

end