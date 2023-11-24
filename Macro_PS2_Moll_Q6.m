% Problem Set 2 (Moll):
% Numerical Dynamic Programming: Growth Model
% by Ed Manuel and Alexander Newton
% This taqkes vfi_IID.m provided on moodle and tweaks it to solve the
% entrepeneur problem in the problem set


%% PARAMETERS

% preferences
risk_aver   = 1;
beta        = 0.95;

%returns
r           = 0.03;
R           = 1.05;

% productivity risk: discretized N(mu,sigma^2)
mu_z    = 3;
sd_z    = 1;
nz      = 5;

% asset grids
na          = 1000;
amax        = 15; 
borrow_lim  = 0.001;
agrid_par   = 1; %1 for linear, 0 for L-shaped

% capital limit
lambda      = 3;

% production function
alpha       = 1/3;
A           = 1;

% computation
max_iter    = 1000;
tol_iter    = 1.0e-6;
Nsim        = 50000;
Tsim        = 500;


%% OPTIONS
Display     = 1;
DoSimulate  = 1;
MakePlots   = 1;

%% DRAW RANDOM NUMBERS
rng(2017);
zrand = rand(Nsim,Tsim);

%% SET UP GRIDS

% assets
agrid = linspace(0,1,na)';
agrid = agrid.^(1./agrid_par);
agrid = borrow_lim + (amax-borrow_lim).*agrid;

% productivity: disretize normal distribution
width = fzero(@(x)discrete_normal(nz,mu_z,sd_z,x),2);
[temp,zgrid,zdist] = discrete_normal(nz,mu_z,sd_z,width);
zcumdist = cumsum(zdist);

%% CONSTRUCT PROFIT GRID

% FOC for k - this is nz by 1 vector
k_uncon=(R./(alpha*A*zgrid)).^(1/(alpha-1));

% Find constrained solution for k and implied profit

% loop over assets
for ia = 1:na
   
    % loop over z
    for iz=1:nz 

   % Constrain k to be less than lambda*a  
   if k_uncon(iz)> lambda*agrid(ia);
   k(ia,iz) = lambda * agrid(ia);
   else 
   k(ia,iz) = k_uncon(iz)  ;  
   end

% Calculate profits given choice of k
Prof(ia,iz) = zgrid(iz) * A * k(ia,iz) ^ alpha - R* k(ia,iz);

        end
    end


%% UTILITY FUNCTION

if risk_aver==1
    u = @(c)log(c);
else    
    u = @(c)(c.^(1-risk_aver)-1)./(1-risk_aver);
end    

%% INITIALIZE VALUE FUNCTION - check if this needs changing

Vguess = zeros(na,nz);
for iz = 1:nz
    Vguess(:,iz) = u(r.*agrid+zgrid(iz))./(1-beta);
end
% Vguess = ones(na,ny);

%% ITERATE ON VALUE FUNCTION

V = Vguess;

Vdiff = 1;
iter = 0;

while iter <= max_iter && Vdiff>tol_iter
    iter = iter + 1;
    Vlast = V;
    V = zeros(na,nz);
    sav = zeros(na,nz);
    savind = zeros(na,nz);
    con = zeros(na,nz);
    
    % loop over assets
    for ia = 1:na
        
        % loop over productivity
        for iz = 1:nz
            cash = Prof(ia,iz) + (1+r)*agrid(ia) ;
            Vchoice = u(max(cash-agrid,1.0e-1000)) + beta.*(Vlast*zdist);          
            [V(ia,iz),savind(ia,iz)] = max(Vchoice);
            sav(ia,iz) = agrid(savind(ia,iz));
            con(ia,iz) = cash - sav(ia,iz);
        end
    end
    
    
    Vdiff = max(max(abs(V-Vlast)));
    if Display >=1
        disp(['Iteration no. ' int2str(iter), ' max val fn diff is ' num2str(Vdiff)]);
    end
end    



%% SIMULATE
if DoSimulate ==1
    zindsim = zeros(Nsim,Tsim);
    aindsim = zeros(Nsim,Tsim);
    
    % initial assets
    aindsim(:,1) = aindsim_final;
    
    %loop over time periods
    for it = 1:Tsim
        if Display >=1 && mod(it,100) ==0
            disp([' Simulating, time period ' int2str(it)]);
        end
        
        %productivity realization: note we vectorize simulations at once because
        %of matlab, in other languages we would loop over individuals
        zindsim(zrand(:,it)<= zcumdist(1),it) = 1;
        for iz = 2:nz
            zindsim(zrand(:,it)> zcumdist(iz-1) & zrand(:,it)<=zcumdist(iz),it) = iz;
        end
        
        % asset choice
        if it<Tsim
            for iz = 1:nz
                aindsim(zindsim(:,it)==iz,it+1) = savind(aindsim(zindsim(:,it)==iz,it),iz);
            end
        end
    end
    
    %assign actual asset and productivity values;
    asim = agrid(aindsim);
    zsim = zgrid(zindsim);


end


%% Consumption dynamics 

% Consumption for each person / point in time
ind = sub2ind(size(con), aindsim, zindsim);
consim=con(ind);

% Average consumption
consim_av= mean(consim,1);

% top-10 share of consumption over time
share_top10 = zeros(size(consim_av))
for t=1:Tsim
    temp=sort(consim(:,t)) ;
    cons_all = sum(temp);
    cons_top10 = sum(temp((Nsim-Nsim*0.1):end));
    share_top10 (t)= cons_top10/cons_all;
end

%% MAKE PLOTS
   
Plot_t=0:50;

figure(1);
    % average consumption over time 
    l1=plot(Plot_t,consim_av(Plot_t+1),'b-','LineWidth',1);
    grid;
%     title('Consumption Policy Function');
    title('Consumption Dyanmics: Average Consumption');

    figure(2);
    %  consumption top 10% over time 
    l1=plot(Plot_t,share_top10(Plot_t+1),'b-','LineWidth',1);
    grid;
%     title('Consumption Policy Function');
    title('Consumption Dyanmics: Share of Top 10');

 
