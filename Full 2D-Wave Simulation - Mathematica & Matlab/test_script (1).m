%% test absolute fitness alone
pop = [1e7-1e5 ; 1e5];
fit_a = [1 0];
fit_r = [1];

t=1:1000;
nt = ones(2,1000);
nt(:,1) = pop;

for i = 2:length(t)
    pop = BM_lottery_model(pop,fit_a,fit_r,T,b,d,sd,sc);
    nt(:,i) = pop;
end

plot(t,nt(1,:),t,nt(2,:))


% test relative fitness alone
pop = [1e7-1e5 ; 1e5];
fit_a = [1]
fit_r = [1]

t=1:1000;
nt = ones(2,1000);
nt(:,1) = pop;

for i = 2:length(t)
    pop = BM_lottery_model(pop,fit_a,fit_r,T,b,d,sd,sc);
    nt(:,i) = pop;
end

plot(t,nt(1,:),t,nt(2,:))

% test relative and absolute fitness

%% test Absolute fitness evolution (no relative fitness) distribution of extinction times versus Te

% basic parameters
T = 1e9;
b = 2;
d = 100/98;
sa = 0.02;
sr = 0.11;
ua = 1e-6;
uad = 0;
ur = 0;
urd = 0;

tEnv = 130;
EnvR = 1./tEnv;
tExt = zeros(size(tEnv));
cutoff = 10/min([sr sa]);

steps = 100000;
collect_data = 1;
outputfile = 'data/compEvo2d_data_ml-test';

% initial population and fitness data
init_pop = [1e7; 1e6];
init_fit_a = [-15 -14];
init_fit_r = [1];

tic
for i = 1:length(tEnv)
    tExt(i) = stochastic_simulation_two_traits_rel_vs_abs( ...
                                    T,init_pop,init_fit_a,init_fit_r,b,d,sa,ua,uad,sr,ur,urd,EnvR(i), ...
                                    steps,collect_data,[outputfile '-' num2str(i)]);
    [i tEnv(i) tExt(i)]
end
toc
dlmwrite('data/compEvo2d_data_ext_times_ml-test.dat',[tEnv' tExt'],'delimiter',',','precision',16);

%% Testing quasi-equilibrium population sizes
T = 1e9;
b = 1;
d = 50/48;
sa = 0.01;
sr = 0.01;
ua = 1e-6;
uad = 1e-5;
ur = 0;
urd = 0;

pop = [1e8];
fit_a = [-42];
fit_r = [1];
endTime = 60;
t=1:endTime;
nt = ones(1,endTime);
nt(:,1) = pop;

for i = 2:endTime
    [pop,delta_ni, Wi,W_bar,wi,c_bar] = BM_lottery_model(nt(:,i-1),fit_a,fit_r,T,b,d,sa,sr);
    nt(:,i) = pop+delta_ni;
end

plot(t,log10(nt(1,:)))

%% Absolute fitness evolution at attractor equilibrium

% basic parameters
% alpha = (0.065/1.85);
T = 1e9; b = 2.0; d = 100/98; sa = 0.01; sr = 0.3;
ua = 1e-6; uad = 1e-5;
% ur = 0;
% urd = 0;
ur = 1e-5;
urd = 1e-5;

rng(5)

tEnv = 130;
EnvR = 1./tEnv;
cutoff = 10/min(sr,sa);

steps = 80000;
collect_data = 1;
start_time = 20000;
end_time = 80000;
outputfile = 'data/compEvo2d_data_ml-222';

% initial population and fitness data
init_pop = [9e8];
init_fit_a = [-1];
init_fit_r = [1];

t_ext = eq_test_stochastic_simulation_two_traits_rel_vs_abs( ...
                                    T,init_pop,init_fit_a,init_fit_r,b,d,sa,ua,uad,sr,ur,urd,EnvR, ...
                                    steps,collect_data,outputfile,start_time,end_time);

%% Plotting selection coefficients for specific examples

% parameters
T = 1e9; b = 2.0; d = 100/98; sa = 0.01; sr = 1.85;
ua = 1e-6; uad = 1e-5; ur = 1e-5; urd = 1e-5;

% looking for a selective advantage of 1.35e-2

% death rates of each class and equilibrium population sizes
d_i = @(k) d./(1-k*sa*d);
N_eq =@(k) 2.*b./(b+2).*T.*(1-(d_i(k)-1)./b);

tend = 1200; t=1:tend; 

% initial population to test selection on relative fitness
pop = [N_eq(8)-10 10];
fit_a = [-48];
fit_r = [1 2];
nt_rel = ones(2,tend);
nt_rel(:,1) = pop;

pt_rel = ones(2,tend);
pt_rel(:,1) = pop/sum(sum(pop));
sel_rel = [];

for i = 2:length(t)
    [newpop,delta_ni] = BM_lottery_model(pop,fit_a,fit_r,T,b,d,sa,sr);
    pop = newpop+delta_ni;
    nt_rel(:,i) = pop';
    pt_rel(:,i) = nt_rel(:,i)/sum(sum(nt_rel(:,i)));
    sel_rel(i-1) = (pt_rel(2,i)-pt_rel(2,i-1))./pt_rel(2,i-1);
end

plot(t(1:end-1),sel_rel)

% initial population parameters for selection on absolute fitnes
pop = [N_eq(31)-10 ; 10];
fit_a = [-31 -30];
fit_r = [1];
nt_abs = ones(2,tend);
nt_abs(:,1) = pop;

pt_abs = ones(2,tend);
pt_abs(:,1) = pop/sum(sum(pop));
sel_abs = [];

for i = 2:length(t)
    [newpop,delta_ni] = BM_lottery_model(pop,fit_a,fit_r,T,b,d,sa,sr);
    pop = newpop+delta_ni;
    nt_abs(:,i) = pop;
    pt_abs(:,i) = nt_abs(:,i)/sum(sum(nt_abs(:,i)));
    sel_abs(i-1) = (pt_abs(2,i)-pt_abs(2,i-1))./pt_abs(2,i-1);
end

plot(t(1:end-1),sel_abs)

plot(t(1:end-1),sel_rel,t(1:end-1),sel_abs)


%% Competition between absolute & relative fitness class

% parameters
T = 1e9; b = 2.0; d = 100/98; sa = 0.01; sr = 0.01;
ua = 1e-6; uad = 1e-5; ur = 1e-5; urd = 1e-5;

% death rates of each class and equilibrium population sizes
d_i = @(k) d./(1-k*sa*d);
N_eq =@(k) 2.*b./(b+2).*T.*(1-(d_i(k)-1)./b);

pop = [N_eq(31)-10 5;5 0];
fit_a = [-31 -30];
fit_r = [1 2];
tend = 2000;

t=1:tend;
nt = ones(2,2,tend);
pt = ones(2,2,tend);
nt(:,:,1) = pop;
pt(:,:,1) = pop/sum(sum(pop));

for i = 2:length(t)
    [newpop,delta_ni] = BM_lottery_model(pop,fit_a,fit_r,T,b,d,sa,sr);
    pop = newpop+delta_ni;
    nt(:,:,i) = pop;
    pt(:,:,i) = nt(:,:,i)/sum(sum(nt(:,:,i)));
end

for i=1:tend
    n11(i) = nt(1,1,i);
    n12(i) = nt(1,2,i);
    n21(i) = nt(2,1,i);
    
    p11(i) = pt(1,1,i);
    p12(i) = pt(1,2,i);
    p21(i) = pt(2,1,i);
end

for i=1:tend-1
    sel_p11(i) = (pt(1,1,i+1)-pt(1,1,i))./pt(1,1,i);
    sel_p12(i) = (pt(1,2,i+1)-pt(1,2,i))./pt(1,2,i);
    sel_p21(i) = (pt(2,1,i+1)-pt(2,1,i))./pt(2,1,i);
end

% figures showing changes in abundances, frequencies and selective strength
% figure(1), plot(t,n11,t,n12,t,n21); legend('n_{11}','n_{12}','n_{21}');
% figure(2), plot(t,p11,t,p12,t,p21); legend('p_{11}','p_{12}','p_{21}');
figure(3), plot(t(1:end-1),sel_p11,t(1:end-1),sel_p12,t(1:end-1),sel_p21); legend('sel-p_{11}','sel-p_{12}','sel-p_{21}'); axis([1 2000 -.1 0.1]);

%% Examining selection coefficients of BM_lotter_model

% parameters
T = 1e9; b = 2.0; d = 100/98; sa = 0.01; sr = 1.85;
ua = 1e-6; uad = 1e-5; ur = 1e-5; urd = 1e-5;

% death rates of each class and equilibrium population sizes
d_i = @(k) d./(1-k*sa*d);
N_eq =@(k) 2.*b./(b+2).*T.*(1-(d_i(k)-1)./b);

pop = [1e8-10 5;5 0];
% pop = [9.6855e-09 5.7984e+08;4.7853e-08 0];
fit_a = [-50 -49];
fit_r = [1 2];
tend = 2000;

t=1:tend;
nt = ones(2,2,tend);
pt = ones(2,2,tend);
nt(:,:,1) = pop;
pt(:,:,1) = pop/sum(sum(pop));

for i = 2:length(t)
    [newpop,delta_ni] = BM_lottery_model(pop,fit_a,fit_r,T,b,d,sa,sr);
    pop = newpop+delta_ni;
    nt(:,:,i) = pop;
    pt(:,:,i) = nt(:,:,i)/sum(sum(nt(:,:,i)));
end

for i=1:tend
    n11(i) = nt(1,1,i);
    n12(i) = nt(1,2,i);
    n21(i) = nt(2,1,i);
    
    p11(i) = pt(1,1,i);
    p12(i) = pt(1,2,i);
    p21(i) = pt(2,1,i);
end

for i=1:tend-1
    sel_p11(i) = (pt(1,1,i+1)-pt(1,1,i))./pt(1,1,i);
    sel_p12(i) = (pt(1,2,i+1)-pt(1,2,i))./pt(1,2,i);
    sel_p21(i) = (pt(2,1,i+1)-pt(2,1,i))./pt(2,1,i);
end

sel_p11(isnan(sel_p11)) = 0;
sel_p12(isnan(sel_p12)) = 0;
sel_p21(isnan(sel_p21)) = 0;

sel_p11(isinf(sel_p11)) = 0;
sel_p12(isinf(sel_p12)) = 0;
sel_p21(isinf(sel_p21)) = 0;

% figures showing changes in abundances, frequencies and selective strength
figure(1), plot(t,n11,t,n12,t,n21); legend('n_{11}','n_{12}','n_{21}');
figure(2), plot(t,p11,t,p12,t,p21); legend('p_{11}','p_{12}','p_{21}');
figure(3), plot(t(1:end-1),sel_p11,t(1:end-1),sel_p12,t(1:end-1),sel_p21); legend('sel-p_{11}','sel-p_{12}','sel-p_{21}'); axis([1 2000 -.1 0.1]);

%% test absolute fitness selection coefficient over absolute fitness space

% this code needs to be fixed!!

% parameters
T = 1e9; b = 2.0; d = 100/98; sa = 0.01; sr = 1.85;
ua = 1e-6; uad = 1e-5; ur = 1e-5; urd = 1e-5;

% death rates of each class and equilibrium population sizes
d_i = @(k) d./(1-k*sa*d);
N_eq =@(k) 2.*b./(b+2).*T.*(1-(d_i(k)-1)./b);

i_ext = floor((1./(sa.*d)).*(1-d./(b+1))); 

i_st = -i_ext+1:-1; 
sela_cff = [];
selr_cff = [];

% compiling selection coefficients for absolute fitness
for i = -i_ext+1:-1
    pop = [N_eq(-i)-10 ; 10];
    fit_a = [-i -i+1];
    fit_r = [1];
    
    for j = 1:200
        [newpop,delta_ni] = BM_lottery_model(pop,fit_a,fit_r,T,b,d,sa,sr);
        if (j<200)
            pop = newpop+delta_ni;
        end
    end
    
    p0 = pop./sum(sum(pop));
    p1 = (newpop+delta_ni)./sum(sum(newpop+delta_ni));
    sela_cff(i+i_ext) = (p1(2)-p0(2))./p0(2);
end

% compiling selection coefficients for relative fitness
for i = -i_ext+1:-1
    pop = [N_eq(-i)-10 10];
    fit_a = [-i];
    fit_r = [1 10];

    for j = 1:200
        [newpop,delta_ni] = BM_lottery_model(pop,fit_a,fit_r,T,b,d,sa,sr);
        if (j<200)
            pop = newpop+delta_ni;
        end
    end
    
    p0 = pop./sum(sum(pop));
    p1 = (newpop+delta_ni)./sum(sum(newpop+delta_ni));
    selr_cff(i+i_ext) = (p1(2)-p0(2))./p0(2);
end

figure(1), plot(i_st,selr_cff,i_st,sela_cff); legend('rel','abs');

%% Checking for Dim.Returns in ci selective advantage

% parameters
T = 1e9; b = 2.0; d = 100/98; sa = 0.01; sr = 1.85;
ua = 1e-6; uad = 1e-5; ur = 1e-5; urd = 1e-5;

% death rates of each class and equilibrium population sizes
d_i = @(k) d.*(1+sa).^-k;
N_eq =@(k) 2.*b./(b+2).*T.*(1-(d_i(k)-1)./b);

i_ext = floor(log((b+1)./d)./log(1+sa)); 
tend = 1200;
t = 1:tend;

% compiling selection coefficients for absolute fitness
for i = 1:1
    pop = [N_eq(-20)-10; 10];
    fit_a = [-20 -18];
    fit_r = [i];
    
    for j = 1:tend
        [newpop,delta_ni] = BM_lottery_model(pop,fit_a,fit_r,T,b,d,sa,sr);
        newpop = newpop+delta_ni;
        newpop
        pi = pop(2)/sum(sum(pop));
        pf = newpop(2)/sum(sum(newpop));
        sela_cff(i,j) = (pf-pi)/pi;
        pop = newpop;
    end
end

figure(1), plot(t(2:end)',sela_cff(1,2:end)'); legend('1','2','3')

%% code for extracting selection coefficients
% parameters
T = 1e9; b = 2.0; d = 100/98; sa = 0.01; sr = 1.85;
ua = 1e-6; uad = 1e-5; ur = 1e-5; urd = 1e-5;

% looking for a selective advantage of 1.35e-2

% death rates of each class and equilibrium population sizes
d_i = @(k) d./(1-k*sa*d);
N_eq =@(k) 2.*b./(b+2).*T.*(1-(d_i(k)-1)./b);

tend = 1200; t=1:tend; 
i_ext = floor((1./(sa.*d)).*(1-d./(b+1))); 

sel_rel = 1:i_ext;
get_sel_rel = true;

Neq = 1:i_ext;
get_Neq = true;

for k = -i_ext:-1
    get_sel_rel = true;
    
    % initial population to test selection on relative fitness
    pop = [N_eq(8)-10 10];
    fit_a = [k];
    fit_r = [1 2];
    nt_rel = ones(2,tend);
    nt_rel(:,1) = pop;

    pt_rel = ones(2,tend);
    pt_rel(:,1) = pop/sum(sum(pop));
    
    s_val = 2;
    
    while((get_sel_rel) && (s_val <= tend))
        
        [newpop,delta_ni] = BM_lottery_model(pop,fit_a,fit_r,T,b,d,sa,sr);
        pop = newpop+delta_ni;
       
        nt_rel(:,s_val) = pop';
        pt_rel(:,s_val) = nt_rel(:,s_val)/sum(sum(nt_rel(:,s_val)));
        i_ext+k+1
        
        if(abs(sel_rel(i_ext+k+1) - (pt_rel(2,s_val)-pt_rel(2,s_val-1))./pt_rel(2,s_val-1))<1e-3)
            sel_rel(i_ext+k+1) = (pt_rel(2,s_val)-pt_rel(2,s_val-1))./pt_rel(2,s_val-1);
            Neq(i_ext+k+1) = sum(sum(pop));
            get_sel_rel = false;
        else
            sel_rel(i_ext+k+1) = (pt_rel(2,s_val)-pt_rel(2,s_val-1))./pt_rel(2,s_val-1);
            Neq(i_ext+k+1) = sum(sum(pop));
            s_val = s_val+1;
        end
    end
end

plot(t(1:end-1),sel_rel)
%%

% parameters (testing parameters)
T = 1e9; b = 2.0; d = 100/98; sa = 0.01; sr = 1.85;
min_diff = 1e-3;
max_time = 1200;

[Neq1,Neq2,Neq3,Neq4,selR1,selR2,selR3,selR4,states]=get_MC_parameters(T,b,d,sa,sr,min_diff,max_time);

figure(1), plot(states,Neq1,states,Neq2,states,Neq3,states,Neq4); legend('numODE','numRoot','apprxOpt','apprxExt','Location','northwest');
figure(2), plot(states,selR1,states,selR2,states,selR3,states,selR4); legend('numODE','numRoot','apprxOpt','apprxExt','Location','south');

