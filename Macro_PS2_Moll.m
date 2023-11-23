% Problem Set 2 (Moll):
% Numerical Dynamic Programming: Growth Model
% by Ed Manuel and Alexander Newton
% This taqkes vfi_IID.m provided on moodle and tweaks it to solve the
% entrepeneur problem in the problem set

clear;
close all;

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
lambda      = 2;

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
Prof(ia,iz) = zgrid(iz) * A * k_con ^ alpha - R* k_con;

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
    aindsim(:,1) = 1;
    
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

%% Find capital and marginal product of capital in final period across individuals
       
% Locate capital for each individual in the k grid

ksim = zeros(Nsim, 1);
       for i = 1:Nsim
   ksim(i,1) = k(aindsim(i,Tsim),zindsim(i,Tsim));
       end


  % Marginal product of capital 
   MPksim = alpha*A.*zsim(:,Tsim).*ksim.^(alpha-1);

%% MAKE PLOTS
if MakePlots ==1 
    figure(1);
        
    % consumption policy function
    subplot(2,4,1);
    plot(agrid,con(:,3),'b-','LineWidth',1);
    grid;
    xlim([0 amax]);
%     title('Consumption Policy Function');
    title('Consumption');

    % savings policy function
    subplot(2,4,2);
    plot(agrid,sav(:,3),'b-',agrid,agrid,'r-','LineWidth',1);
    hold on;
    plot(agrid,zeros(na,1),'k','LineWidth',0.5);
    hold off;
    grid;
    xlim([0 amax]);
%     title('Savings Policy Function (a''-a)');
    title('Savings');
    
    % consumption policy function: zoomed in
    % subplot(2,4,3);
    % plot(agrid,con(:,3),'b-','LineWidth',2);
    % grid;
    % xlim([0 5]);
    % title('Consumption: Zoomed');
    
     % savings policy function: zoomed in
    subplot(2,4,3:4);
    plot(agrid,sav(:,3),'b-',agrid,agrid,'r-','LineWidth',2);
    hold on;
    plot(agrid,zeros(na,1),'k','LineWidth',0.5);
    hold off;
    grid;
    xlim([0 5]);
    title('Savings: Zoomed ');
    
    
    %Prod distribution
    % subplot(2,4,5);
    % hist(zsim(:,Tsim),zgrid);
    % h = findobj(gca,'Type','patch');
    % set(h,'FaceColor',[0 0.5 0.5],'EdgeColor','blue','LineStyle','-');
    % ylabel('')
    % title('Productivity distribution');
    % 
    %asset distribution
    subplot(2,4,5:6);
    hist(asim(:,Tsim),[0:0.05:10]);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',[.7 .7 .7],'EdgeColor','black','LineStyle','-');
    ylabel('')
    title('Asset distribution');

        %Marginal product of cqpital distribution
    subplot(2,4,7:8);
    hist(MPksim,[0:0.001:3]);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',[.7 .7 .7],'EdgeColor','black','LineStyle','-');
    ylabel('')
    title('Marginal Product of Capital distribution');

    %convergence check
    % subplot(2,4,8);
    % plot([1:Tsim]',mean(asim,1),'k-','LineWidth',1.5);
    % ylabel('Time Period');
    % title('Mean Asset Convergence');
    % 
    
   % asset distribution statistics
    aysim = asim(:,Tsim) ./ mean(zsim(:,Tsim));
    disp(['Mean assets: ' num2str(mean(aysim))]);
    disp(['Fraction borrowing constrained: ' num2str(sum(aysim==borrow_lim)./Nsim * 100) '%']);
    disp(['10th Percentile: ' num2str(quantile(aysim,.1))]);
    disp(['50th Percentile: ' num2str(quantile(aysim,.5))]);
    disp(['90th Percentile: ' num2str(quantile(aysim,.9))]);
    disp(['99th Percentile: ' num2str(quantile(aysim,.99))]);

end
