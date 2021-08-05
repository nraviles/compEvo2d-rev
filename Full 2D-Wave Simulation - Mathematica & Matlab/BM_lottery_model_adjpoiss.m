function [newpop,delta_ni, Wi,W_bar,wi,c_bar] = BM_lottery_model(pop,fit_a,fit_r,T,b,d,sa,sr)
% BM_lottery_model takes in the set of population abundances for each class
% and the fitness array to calculate the the new expected abundances. The
% model is specified in Bertram and Masel 2019 - "Density-dependent 
% selection and the limits of relative fitness" published in TPB. Mutations
% are also incorporated in the function since mutations are determined by
% the number of births.
%
% inputs:
% pop = array with current set of abundances 
% fit_a = number of mutations in absolute fitness trait (1st index in pop)
% fit_r = number of mutations in relative fitness trait (2nd index in pop)
% %
% T = total available territories
% b = base birth rate 
% d = base death rate 
% sd = selection coefficient for beneficial mutation descreasing death rate
% sr = selection coefficient for beneficial mutation increasing competitive ability 
% 
% 
% outputs:
% new_pop = expected abundances due to selection
% delta_ni = change to class ni
% Wi = set of absolute fitness values
% W_bar = mean absolute fitness
% wi = set of relative fitness values (mean of wi = 1)
%

%digits(16); % I commented this out

% ------------------- Selection ------------------------------------------

% basic dimensions
ka = size(fit_a,2);
kr = size(fit_r,2);

% To avoid numerical error involved with abundances < 1e-4, I set any 
% abundance smaller than 1e-3 to be 0, since the probability that such a
% class remains is small according to a poisson distribution used in the
% main model.
pop(pop<1e-2) = 0;

% expected change \Delta_+ n_i
Na = sum(sum(pop));     % current population size
U = T - Na;             % unoccupied territories 

mi= pop.*(b.*U./T);     % array of avg propagules dispersed per class
li= pop.*(b./T);

L = sum(sum(li));
ci = ones(ka,1)*((1+sr).^fit_r);       % relative fitness is (1 + sr)^#mut x
c_bar = sum(sum(mi.*ci))./(sum(sum(mi)));

% IMPORTANT!! the ci should never be zero or negative, since this would lead
% to problems in calculating the term (Ri + Ai).*ci./c_bar below for
% births. This could happen if there are deleterious mutations in the "c"
% trait --- I include in the mutation function
% two_dim_abs_vs_rel_mutations function --- and "sr > 1."

inv_di = ((1./(d.*(1+sa).^(-fit_a)))').*ones(1,kr);  % calculate inverse because can't have di-->inf
inv_di(inv_di<0) = 0;                           % inverse of death rate should not be negative, set these values to zero 

% NOTE: with inv_di defined above, 1/d_i-1 -1/di = sa, 0<= 1/di <= 1/d, 1 <= di < inf;

% Technically this birthing process is improper, in terms of straight births before competition and such is just a random sample of Poiss(b * n_ij) for each pop
% births = poissrnd(b .* pop) I believe from here delta_ni is not just the ind_di .* pop and ind_di .* births but rather the result of competition
% let us try and conver Jasons code for competitive wins to matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% New Jason %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def deltnplussim(m,c,U): % here we will need m =  births or maybe better said as propogules, and here U is just a number corresponding to T - N, the unoccupied territories
#    scatter=np.zeros([int(U),len(m)])
#    for i in range(len(m)):
#        for y in xrange(int(m[i])):
#            cell=int(int(U)*np.random.rand());
#            scatter[cell,i]=scatter[cell,i]+1;

    l=m/float(U)        
    scatter=np.random.poisson(lam=l,size=[U,len(m)]) #dispersion of propugules 
    
    winsnocomp=np.zeros(len(m)); wins1=np.zeros(len(m)); wins2=np.zeros(len(m)); %here it gets weird, he attempts to account for the conditional probabilities, but in true all
	% we need here is who wins. this will be odd as m, the propagules is now an array (i,j) from here its as simple as finding the index of the propagules that scattered to the
	% same location and then using cumsum and a random sampling to identify a winner. From here I believe we do not need to calculate the types of wins, but rather just WHO won
	% I have no clue where to start with this.
    comp=np.zeros(U);
    for i in range(int(U)):
        comp[i]=sum(scatter[i]) #total number competing per territory
        if comp[i]>0:            
            lotterycmf=np.cumsum(np.array(scatter[i])*c) # Sum mi ci / Sum mi, cbar n *c, n*c + m*c, ..., Sum mi ci is lotterycmf[-1]
            victor=bisect.bisect(lotterycmf,np.random.rand()*lotterycmf[-1]) #random.rand random between 0-1, [0 , c1 m1, c1 m1 + c2 m2] winner based on uniform
            
            if scatter[i][victor]==1 and comp[i]==1:
                winsnocomp[victor]=winsnocomp[victor]+1
            elif scatter[i][victor]==1:
                wins1[victor]=wins1[victor]+1
            else: 
                wins2[victor]=wins2[victor]+1
        # wins based on size: wins based on no comp, based on 1-on-1, 1-on-many
    return np.array([winsnocomp,wins1,wins2])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% New Jason %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
newpop = (inv_di).*pop;             % left over adults after death
delta_ni = (inv_di).*births;        % juviniles that make it to be adults
% NOTE!!! the total population is newpop+delta_ni, but here I keep them
% seperate becuase only the delta_ni portion experiences mutations.

Wi = (newpop + delta_ni)./pop;
W_bar = sum(sum(newpop + delta_ni))/Na;
wi = Wi./W_bar;

end