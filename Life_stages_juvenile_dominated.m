%{
This code has no pre-allocation.
Completely deterministic population dynamics, only stochastic part is
drawing of mutants from populations.
Population densities check precisely with the ones in Mathematica
Last changed: 2021-07-05 -> corrected pfix algorithm
%}

function Life_stages_juvenile_dominated

clf
%clc
clear
rng('shuffle')

% Parameters
time    = 5e2; % generations
dt      = 1e-3; % ODE integration step size
lim     = 0.01; % pop extinction

za        = -1.5; % trade-off
zj        = -1; % trade-off
Tgen      = 0.001; % generalist phenotype
jmax      = (1+(-1).*exp(1).^((-0.5E0).*zj)).^(-1).*(Tgen+(-1).*exp(1).^((-1) .*zj).*Tgen);
amax      = (1+(-1).*exp(1).^((-0.5E0).*za)).^(-1).*(Tgen+(-1).*exp(1).^((-1) .*za).*Tgen);

aA1       = 0.5;
aA2       = 0.5;
aJ1       = 0.25;
aJ2       = 0.25;
dA        = 1;
dJ        = 0.6;
KA        = 3000;
KJ        = 1500;
rA1        = 6;
rA2        = 6;
rJ1        = 6;
rJ2        = 6;

sc        = 0.02; % variance
sp        = 0.015; % fixed step size
xJ_state  = [];
xA_state  = [];
CJ_state  = [];
CA_state  = [];
eJ1_state = [];
eJ2_state = [];
eA1_state = [];
eA2_state = [];

sk = nan;
parent = nan;
mu = nan;

% Saving data
pname = '/Users/prm/Documents/MATLAB_scripts/Project2/Deterministic/isolated_mutations/Data';
dname = ['JD_symmetric_zj_' num2str(zj) '_za' num2str(za) '_' date];
mkdir(fullfile(pname,dname));
Data = {};

% Create data structure for saving data
create_data

% Plotting specs
colres  = 64;
col     = [max(min(3/2-abs(3-(1:colres)*4/colres),1),0)' max(min(3/2-abs(2-(1:colres)*4/colres),1),0)' max(min(3/2-abs(1-(1:colres)*4/colres),1),0)'];
figure(1);
clf(1);
set(gcf,'Position',[100 100 700 700])

% Simulation
for iter = 1:2
    iter
    clf
    
    % Initial conditions
    RJ1     = KJ;
    RJ2     = KJ;
    RA1     = KA;
    RA2     = KA;
    CJ      = 100;
    CA      = 100;
    xA      = 0.4;
    xJ      = 0.4;
    xA_in   = xA;
    xJ_in   = xJ;
    xJ_parent = nan;
    xA_parent = nan;
    xJm1    = nan;
    xJm2    = nan;
    xAm1    = nan;
    xAm2    = nan;
    
    eA1     = amax * ((1 - exp(1) .^ (-za .* (1 - xA))) ./ (1 - exp(1) .^ (-za)));
    eA2     = amax * ((1 - exp(1) .^ (-za .*      xA))  ./ (1 - exp(1) .^ (-za)));
    eJ1     = jmax * ((1 - exp(1) .^ (-zj .* (1 - xJ))) ./ (1 - exp(1) .^ (-zj)));
    eJ2     = jmax * ((1 - exp(1) .^ (-zj .*      xJ))  ./ (1 - exp(1) .^ (-zj)));
    
    filename = fullfile(pname,dname,['iter_' num2str(iter) '_data.mat'])
    filename_fig = fullfile(pname,dname,['iter_' num2str(iter) '.pdf'])
    
    
    for t   = 1:time
        %% Population dynamics of resources and consumers
        euler
        
        %% Extinction of dead cells
        extinction
        
        %% Save data
        save_data
        
        %% Plotting
        %plotting
        %fprintf('time %2.0f, r = %2.0f \n',t,length(CA))
        
        %% Mutations
        choose_parent
        variable_step_mutation
        switch_mu
        
        %% Merge populations with repeated phenotypes
        merge_pops
        
        % Test for time out
        if sk == 1
            sj = 1;
            continue
        else
            sj = 0;
        end
    end
    
    if sj == 1
        %% Saving figure
        %saveas(gcf,filename_fig);
        continue
    end
    
    %% Saving figure
    %saveas(gcf,filename_fig);
end

%% NESTED FUNCTIONS
    function create_data
        Data   = cell(time+1,4);
        Data{1,1} = "Juveniles";
        Data{1,2} = "Juvenile phenotypes";
        Data{1,3} = "Adults";
        Data{1,4} = "Adult phenotypes";
    end

    function euler
        dRJ1_dt = inf;
        dRJ2_dt = inf;
        dRA1_dt = inf;
        dRA2_dt = inf;
        dCJ_dt = inf;
        dCA_dt = inf;
        max_abs_change = 1e-6; % maximum change size for Euler's forward method
        tic
        while abs(dRJ1_dt) > max_abs_change | abs(dRJ2_dt) > max_abs_change | abs(dRA1_dt) > max_abs_change | abs(dRA2_dt) > max_abs_change | abs(dCJ_dt) > max_abs_change | abs(dCA_dt) > max_abs_change
            dRJ1_dt = (rJ1 * RJ1 * (1 - RJ1 / KJ) - sum(eJ1 .* RJ1 .*CJ));
            RJ1    = RJ1 +  dRJ1_dt * dt;
            dRJ2_dt = (rJ2 * RJ2 * (1 - RJ2 / KJ) - sum(eJ2 .* RJ2 .*CJ));
            RJ2    = RJ2 +  dRJ2_dt * dt;
            
            dRA1_dt = (rA1 * RA1 * (1 - RA1 / KA) - sum(eA1 .* RA1 .*CA));
            RA1    = RA1 +  dRA1_dt * dt;
            dRA2_dt = (rA2 * RA2 * (1 - RA2 / KA) - sum(eA2 .* RA2 .*CA));
            RA2    = RA2 +  dRA2_dt * dt;
            
            dCJ_dt = CA .* (aA1 .* eA1 .* RA1 + aA2 .* eA2 .* RA2)+(-1) .* CJ .* (dJ + aJ1 .* eJ1 .* RJ1+ ...
                aJ2 .* eJ2 .* RJ2);
            CJ    = CJ +  dCJ_dt * dt;
            dCA_dt = (-1) .* CA .* dA + aJ1 .* CJ .* eJ1 .* RJ1+aJ2 .* CJ .* eJ2 .* RJ2;
            CA    = CA +  dCA_dt * dt;
            
            time_cond = toc;
            if time_cond >= 180
                sk = 1;
                break
            else
                sk = 0;
            end
        end
    end

    function extinction
        if any(CJ <= lim) || any(CA <= lim)
            ej      = find(CJ <= lim);
            CJ(ej)  = [];
            CA(ej)  = [];
            xJ(ej)  = [];
            xA(ej)  = [];
            eJ1(ej) = [];
            eJ2(ej) = [];
            eA1(ej) = [];
            eA2(ej) = [];
            ea      = find(CA <= lim);
            CJ(ea)  = [];
            CA(ea)  = [];
            xJ(ea)  = [];
            xA(ea)  = [];
            eJ1(ea) = [];
            eJ2(ea) = [];
            eA1(ea) = [];
            eA2(ea) = [];
        end
    end

    function plotting
        figure(1)
        scatter(xJ, xA, sqrt(CJ/sum(CJ))*2e2, 'MarkerEdgeColor', col(ceil(t/time*colres),:))
        hold on;
        xlabel('x_J')
        ylabel('x_A')
        axis([0 1 0 1])
        title({['Juvenile dominated. Parameters:']
            ['zj = ', num2str(zj), ', za = ',num2str(za), ', initial \Theta_J = ', ...
            num2str(xJ_in), ', initial \Theta_A = ', num2str(xA_in), ', generations = ', num2str(time) %, ', rJ = ', num2str(rJ), ', KJ = ' num2str(KJ)
            ]})
        
        drawnow;
    end


    function choose_parent
        p = zeros(20,1);
        p(1:numel(CA)) = CA .* (aA1 .* eA1 .* RA1 + aA2 .* eA2 .* RA2);
        p0     = sum(p);
        p1 = sum(rand*p0 <= cumsum(p));
        switch p1
            case 20
                parent = 1;
            case 19
                parent = 2;
            case 18
                parent = 3;
            case 17
                parent = 4;
            case 16
                parent = 5;
            case 15
                parent = 6;
            case 14
                parent = 7;
            case 13
                parent = 8;
            case 12
                parent = 9;
            case 11
                parent = 10;
            case 10
                parent = 11;
            case 9
                parent = 12;
            case 8
                parent = 13;
            case 7
                parent = 14;
            case 6
                parent = 15;
            case 5
                parent = 16;
            case 4
                parent = 17;
            case 3
                parent = 18;
            case 2
                parent = 19;
            case 1
                parent = 20;
        end
    end

    function variable_step_mutation
        xJ_parent = xJ(parent);
        xJm = normrnd(xJ_parent,sc,20,1); % draw 10 mutations with sc as variance and xJ_parent as mean
        xJm1 = xJm(find(xJm > xJ_parent,1)); % pick one that is larger than xJ_parent
        xJm2 = xJm(find(xJm < xJ_parent,1)); % pick one that is smaller than xJ_parent
        xJm1 = min(max(xJm1,0),1); % bind mutation between 0 and 1
        xJm2 = min(max(xJm2,0),1); % bind mutation between 0 and 1
        
        xA_parent = xA(parent);
        xAm = normrnd(xA_parent,sc,20,1); % draw mutation with sc as variance and xA_parent as mean
        xAm1 = xAm(find(xAm > xA_parent,1)); % pick one that is larger than xA_parent
        xAm2 = xAm(find(xAm < xA_parent,1)); % pick one that is smaller than xA_parent
        xAm1 = min(max(xAm1,0),1); % bind mutation between 0 and 1
        xAm2 = min(max(xAm2,0),1); % bind mutation between 0 and 1
        
        eJ1m1   = jmax * ((1 - exp(1) .^ (-zj .* (1 - xJm1))) ./ (1 - exp(1) .^ (-zj)));
        eJ2m1   = jmax * ((1 - exp(1) .^ (-zj .*      xJm1))  ./ (1 - exp(1) .^ (-zj)));
        eJ1m2   = jmax * ((1 - exp(1) .^ (-zj .* (1 - xJm2))) ./ (1 - exp(1) .^ (-zj)));
        eJ2m2   = jmax * ((1 - exp(1) .^ (-zj .*      xJm2))  ./ (1 - exp(1) .^ (-zj)));
        eA1m1   = amax * ((1 - exp(1) .^ (-za .* (1 - xAm1))) ./ (1 - exp(1) .^ (-za)));
        eA2m1   = amax * ((1 - exp(1) .^ (-za .*      xAm1))  ./ (1 - exp(1) .^ (-za)));
        eA1m2   = amax * ((1 - exp(1) .^ (-za .* (1 - xAm2))) ./ (1 - exp(1) .^ (-za)));
        eA2m2   = amax * ((1 - exp(1) .^ (-za .*      xAm2))  ./ (1 - exp(1) .^ (-za)));
        
        b        = (aA1 .* eA1(parent) .* RA1 + aA2 .* eA2(parent) .* RA2)./((aA1 .* eA1(parent) .* RA1 + aA2 .* eA2(parent) .* RA2) + dA); % reproduction prob. res.
        m        = (aJ1 .* eJ1(parent) .* RJ1 + aJ2 .* eJ2(parent) .* RJ2)./((aJ1 .* eJ1(parent) .* RJ1 + aJ2 .* eJ2(parent) .* RJ2) + dJ); %  maturation prob. res.
        b_m1      = (aA1 .* eA1m1 .* RA1 + aA2 .* eA2m1 .* RA2)./((aA1 .* eA1m1 .* RA1 + aA2 .* eA2m1 .* RA2) + dA); % reproduction prob. mut.
        m_m1      = (aJ1 .* eJ1m1 .* RJ1 + aJ2 .* eJ2m1 .* RJ2)./((aJ1 .* eJ1m1 .* RJ1 + aJ2 .* eJ2m1 .* RJ2) + dJ); % maturation prob. mut.
        b_m2      = (aA1 .* eA1m2 .* RA1 + aA2 .* eA2m2 .* RA2)./((aA1 .* eA1m2 .* RA1 + aA2 .* eA2m2 .* RA2) + dA); % reproduction prob. mut.
        m_m2      = (aJ1 .* eJ1m2 .* RJ1 + aJ2 .* eJ2m2 .* RJ2)./((aJ1 .* eJ1m2 .* RJ1 + aJ2 .* eJ2m2 .* RJ2) + dJ); % maturation prob. mut.
        juv_m1 = floor(b .* (1 + m_m1));
        ad_m1  = floor(b_m1 .* (1 + m));
        juv_m2 = floor(b .* (1 + m_m2));
        ad_m2  = floor(b_m2 .* (1 + m));
        pfix_j1 = (m_m1 - ((1-b)/b)).*juv_m1;
        pfix_a1 = (m - ((1 - b_m1)/b_m1)).*ad_m1;
        pfix_j2 = (m_m2 - ((1-b)/b)).*juv_m2;
        pfix_a2 = (m - ((1 - b_m2)/b_m2)).*ad_m2;
        
        probs = [pfix_j1, pfix_a1, pfix_j2, pfix_a2];
        probs_0 = sum(probs);
        mu = sum(rand*probs_0 <= cumsum(probs));
    end

    function switch_mu
        switch mu
            case 4
                if xJm1 ~= 1 || xJm1 ~= 0
                    eJ1m   = jmax * ((1 - exp(1) .^ (-zj .* (1 - xJm1))) ./ ...
                        (1 - exp(1) .^ (-zj))); % Updating trait eJ1
                    eJ2m   = jmax * ((1 - exp(1) .^ (-zj .*      xJm1))  ./ ...
                        (1 - exp(1) .^ (-zj))); % Updating trait eJ2
                    eA1m   = amax * ((1 - exp(1) .^ (-za .* (1 - xA_parent))) ./ ...
                        (1 - exp(1) .^ (-za))); % Updating trait eA1
                    eA2m   = amax * ((1 - exp(1) .^ (-za .*      xA_parent))  ./ ...
                        (1 - exp(1) .^ (-za))); % Updating trait eA2
                    CJm = [1];
                    CAm = [1];
                    
                    xJ  = [xJ; xJm1];
                    xA  = [xA; xA_parent];
                    eJ1 = [eJ1; eJ1m];
                    eJ2 = [eJ2; eJ2m];
                    eA1 = [eA1; eA1m];
                    eA2 = [eA2; eA2m];
                    CJ  = [CJ; CJm];
                    CA  = [CA; CAm];
                else
                end
            case 3
                if xAm1 ~= 1 || xAm1 ~= 0
                    eJ1m   = jmax * ((1 - exp(1) .^ (-zj .* (1 - xJ_parent))) ./ ...
                        (1 - exp(1) .^ (-zj))); % Updating trait eJ1
                    eJ2m   = jmax * ((1 - exp(1) .^ (-zj .*      xJ_parent))  ./ ...
                        (1 - exp(1) .^ (-zj))); % Updating trait eJ2
                    eA1m   = amax * ((1 - exp(1) .^ (-za .* (1 - xAm1))) ./ ...
                        (1 - exp(1) .^ (-za))); % Updating trait eA1
                    eA2m   = amax * ((1 - exp(1) .^ (-za .*      xAm1))  ./ ...
                        (1 - exp(1) .^ (-za))); % Updating trait eA2
                    CJm = [1];
                    CAm = [1];
                    
                    xJ  = [xJ; xJ_parent];
                    xA  = [xA; xAm1];
                    eJ1 = [eJ1; eJ1m];
                    eJ2 = [eJ2; eJ2m];
                    eA1 = [eA1; eA1m];
                    eA2 = [eA2; eA2m];
                    CJ  = [CJ; CJm];
                    CA  = [CA; CAm];
                else
                end
            case 2
                if xJm2 ~= 1 || xJm2 ~= 0
                    eJ1m   = jmax * ((1 - exp(1) .^ (-zj .* (1 - xJm2))) ./ ...
                        (1 - exp(1) .^ (-zj))); % Updating trait eJ1
                    eJ2m   = jmax * ((1 - exp(1) .^ (-zj .*      xJm2))  ./ ...
                        (1 - exp(1) .^ (-zj))); % Updating trait eJ2
                    eA1m   = amax * ((1 - exp(1) .^ (-za .* (1 - xA_parent))) ./ ...
                        (1 - exp(1) .^ (-za))); % Updating trait eA1
                    eA2m   = amax * ((1 - exp(1) .^ (-za .*      xA_parent))  ./ ...
                        (1 - exp(1) .^ (-za))); % Updating trait eA2
                    CJm = [1];
                    CAm = [1];
                    
                    xJ  = [xJ; xJm2];
                    xA  = [xA; xA_parent];
                    eJ1 = [eJ1; eJ1m];
                    eJ2 = [eJ2; eJ2m];
                    eA1 = [eA1; eA1m];
                    eA2 = [eA2; eA2m];
                    CJ  = [CJ; CJm];
                    CA  = [CA; CAm];
                else
                end
            case 1
                if xAm2 ~= 1 || xAm2 ~= 0
                    eJ1m   = jmax * ((1 - exp(1) .^ (-zj .* (1 - xJ_parent))) ./ ...
                        (1 - exp(1) .^ (-zj))); % Updating trait eJ1
                    eJ2m   = jmax * ((1 - exp(1) .^ (-zj .*      xJ_parent))  ./ ...
                        (1 - exp(1) .^ (-zj))); % Updating trait eJ2
                    eA1m   = amax * ((1 - exp(1) .^ (-za .* (1 - xAm2))) ./ ...
                        (1 - exp(1) .^ (-za))); % Updating trait eA1
                    eA2m   = amax * ((1 - exp(1) .^ (-za .*      xAm2))  ./ ...
                        (1 - exp(1) .^ (-za))); % Updating trait eA2
                    CJm = [1];
                    CAm = [1];
                    
                    xJ  = [xJ; xJ_parent];
                    xA  = [xA; xAm2];
                    eJ1 = [eJ1; eJ1m];
                    eJ2 = [eJ2; eJ2m];
                    eA1 = [eA1; eA1m];
                    eA2 = [eA2; eA2m];
                    CJ  = [CJ; CJm];
                    CA  = [CA; CAm];
                else
                end
        end
    end

    function variable_binned_step_mutation
        xJ_parent = xJ(parent);
        xJm = normrnd(xJ_parent,sc,20,1); % draw 10 mutations with sc as variance and xJ_parent as mean
        xJm1 = xJm(find(xJm > xJ_parent,1)); % pick one that is larger than xJ_parent
        xJm2 = xJm(find(xJm < xJ_parent,1)); % pick one that is smaller than xJ_parent
        xJm1 = floor(xJm1/sc)*sc; % locate mutation within one of 50 bins and sets the mutation as that bin's value
        xJm1 = xJm1+0.5*sc; % center mutation around parent phenotype
        xJm1 = min(max(xJm1,0),1); % bind mutation between 0 and 1
        xJm2 = floor(xJm2/sc)*sc; % locate mutation within one of 50 bins and sets the mutation as that bin's value
        xJm2 = xJm2+0.5*sc; % center mutation around parent phenotype
        xJm2 = min(max(xJm2,0),1); % bind mutation between 0 and 1
        
        xA_parent = xA(parent);
        xAm = normrnd(xA_parent,sc,20,1); % draw mutation with sc as variance and xA_parent as mean
        xAm1 = xAm(find(xAm > xA_parent,1)); % pick one that is larger than xA_parent
        xAm2 = xAm(find(xAm < xA_parent,1)); % pick one that is smaller than xA_parent
        xAm1 = floor(xAm1/sc)*sc; % locate mutation within one of 50 bins and sets the mutation as that bin's value
        xAm1 = xAm1+0.5*sc; % center mutation around parent phenotype
        xAm1 = min(max(xAm1,0),1); % bind mutation between 0 and 1
        xAm2 = floor(xAm2/sc)*sc; % locate mutation within one of 50 bins and sets the mutation as that bin's value
        xAm2 = xAm2+0.5*sc; % center mutation around parent phenotype
        xAm2 = min(max(xAm2,0),1); % bind mutation between 0 and 1
        
        eJ1m1   = jmax * ((1 - exp(1) .^ (-zj .* (1 - xJm1))) ./ (1 - exp(1) .^ (-zj)));
        eJ2m1   = jmax * ((1 - exp(1) .^ (-zj .*      xJm1))  ./ (1 - exp(1) .^ (-zj)));
        eJ1m2   = jmax * ((1 - exp(1) .^ (-zj .* (1 - xJm2))) ./ (1 - exp(1) .^ (-zj)));
        eJ2m2   = jmax * ((1 - exp(1) .^ (-zj .*      xJm2))  ./ (1 - exp(1) .^ (-zj)));
        eA1m1   = amax * ((1 - exp(1) .^ (-za .* (1 - xAm1))) ./ (1 - exp(1) .^ (-za)));
        eA2m1   = amax * ((1 - exp(1) .^ (-za .*      xAm1))  ./ (1 - exp(1) .^ (-za)));
        eA1m2   = amax * ((1 - exp(1) .^ (-za .* (1 - xAm2))) ./ (1 - exp(1) .^ (-za)));
        eA2m2   = amax * ((1 - exp(1) .^ (-za .*      xAm2))  ./ (1 - exp(1) .^ (-za)));
        
        b        = (aA1 .* eA1(parent) .* RA1 + aA2 .* eA2(parent) .* RA2)./((aA1 .* eA1(parent) .* RA1 + aA2 .* eA2(parent) .* RA2) + dA); % reproduction prob. res.
        m        = (aJ1 .* eJ1(parent) .* RJ1 + aJ2 .* eJ2(parent) .* RJ2)./((aJ1 .* eJ1(parent) .* RJ1 + aJ2 .* eJ2(parent) .* RJ2) + dJ); %  maturation prob. res.
        b_m1      = (aA1 .* eA1m1 .* RA1 + aA2 .* eA2m1 .* RA2)./((aA1 .* eA1m1 .* RA1 + aA2 .* eA2m1 .* RA2) + dA); % reproduction prob. mut.
        m_m1      = (aJ1 .* eJ1m1 .* RJ1 + aJ2 .* eJ2m1 .* RJ2)./((aJ1 .* eJ1m1 .* RJ1 + aJ2 .* eJ2m1 .* RJ2) + dJ); % maturation prob. mut.
        b_m2      = (aA1 .* eA1m2 .* RA1 + aA2 .* eA2m2 .* RA2)./((aA1 .* eA1m2 .* RA1 + aA2 .* eA2m2 .* RA2) + dA); % reproduction prob. mut.
        m_m2      = (aJ1 .* eJ1m2 .* RJ1 + aJ2 .* eJ2m2 .* RJ2)./((aJ1 .* eJ1m2 .* RJ1 + aJ2 .* eJ2m2 .* RJ2) + dJ); % maturation prob. mut.
        juv_m1 = floor(b .* (1 + m_m1));
        ad_m1  = floor(b_m1 .* (1 + m));
        juv_m2 = floor(b .* (1 + m_m2));
        ad_m2  = floor(b_m2 .* (1 + m));
        pfix_j1 = (m_m1 - ((1-b)/b)).*juv_m1;
        pfix_a1 = (m - ((1 - b_m1)/b_m1)).*ad_m1;
        pfix_j2 = (m_m2 - ((1-b)/b)).*juv_m2;
        pfix_a2 = (m - ((1 - b_m2)/b_m2)).*ad_m2;
        
        probs = [pfix_j1, pfix_a1, pfix_j2, pfix_a2];
        probs_0 = sum(probs);
        mu = sum(rand*probs_0 <= cumsum(probs));
    end

    function fixed_step_mutation
        xJ_parent = xJ(parent);
        xJm1 = xJ_parent + sp;
        xJm2 = xJ_parent - sp;
        xJm1 = min(max(xJm1,0),1); % binds mutation between 0 and 1
        xJm2 = min(max(xJm2,0),1); % binds mutation between 0 and 1
        
        xA_parent  = xA(parent);
        xAm1 = xA_parent + sp;
        xAm2 = xA_parent - sp;
        xAm1 = min(max(xAm1,0),1); % binds mutation between 0 and 1
        xAm2 = min(max(xAm2,0),1); % binds mutation between 0 and 1
        
        eJ1m1   = jmax * ((1 - exp(1) .^ (-zj .* (1 - xJm1))) ./ (1 - exp(1) .^ (-zj)));
        eJ2m1   = jmax * ((1 - exp(1) .^ (-zj .*      xJm1))  ./ (1 - exp(1) .^ (-zj)));
        eJ1m2   = jmax * ((1 - exp(1) .^ (-zj .* (1 - xJm2))) ./ (1 - exp(1) .^ (-zj)));
        eJ2m2   = jmax * ((1 - exp(1) .^ (-zj .*      xJm2))  ./ (1 - exp(1) .^ (-zj)));
        eA1m1   = amax * ((1 - exp(1) .^ (-za .* (1 - xAm1))) ./ (1 - exp(1) .^ (-za)));
        eA2m1   = amax * ((1 - exp(1) .^ (-za .*      xAm1))  ./ (1 - exp(1) .^ (-za)));
        eA1m2   = amax * ((1 - exp(1) .^ (-za .* (1 - xAm2))) ./ (1 - exp(1) .^ (-za)));
        eA2m2   = amax * ((1 - exp(1) .^ (-za .*      xAm2))  ./ (1 - exp(1) .^ (-za)));
        
        b        = (aA1 .* eA1(parent) .* RA1 + aA2 .* eA2(parent) .* RA2)./((aA1 .* eA1(parent) .* RA1 + aA2 .* eA2(parent) .* RA2) + dA); % reproduction prob. res.
        m        = (aJ1 .* eJ1(parent) .* RJ1 + aJ2 .* eJ2(parent) .* RJ2)./((aJ1 .* eJ1(parent) .* RJ1 + aJ2 .* eJ2(parent) .* RJ2) + dJ); %  maturation prob. res.
        b_m1      = (aA1 .* eA1m1 .* RA1 + aA2 .* eA2m1 .* RA2)./((aA1 .* eA1m1 .* RA1 + aA2 .* eA2m1 .* RA2) + dA); % reproduction prob. mut.
        m_m1      = (aJ1 .* eJ1m1 .* RJ1 + aJ2 .* eJ2m1 .* RJ2)./((aJ1 .* eJ1m1 .* RJ1 + aJ2 .* eJ2m1 .* RJ2) + dJ); % maturation prob. mut.
        b_m2      = (aA1 .* eA1m2 .* RA1 + aA2 .* eA2m2 .* RA2)./((aA1 .* eA1m2 .* RA1 + aA2 .* eA2m2 .* RA2) + dA); % reproduction prob. mut.
        m_m2      = (aJ1 .* eJ1m2 .* RJ1 + aJ2 .* eJ2m2 .* RJ2)./((aJ1 .* eJ1m2 .* RJ1 + aJ2 .* eJ2m2 .* RJ2) + dJ); % maturation prob. mut.
        juv_m1 = floor(b .* (1 + m_m1));
        ad_m1  = floor(b_m1 .* (1 + m));
        juv_m2 = floor(b .* (1 + m_m2));
        ad_m2  = floor(b_m2 .* (1 + m));
        pfix_j1 = (m_m1 - ((1-b)/b)).*juv_m1;
        pfix_a1 = (m - ((1 - b_m1)/b_m1)).*ad_m1;
        pfix_j2 = (m_m2 - ((1-b)/b)).*juv_m2;
        pfix_a2 = (m - ((1 - b_m2)/b_m2)).*ad_m2;
        
        probs = [pfix_j1, pfix_a1, pfix_j2, pfix_a2];
        probs_0 = sum(probs);
        mu = sum(rand*probs_0 <= cumsum(probs));
    end

    function save_data
        Data{t+1,1} = CJ;
        Data{t+1,2} = xJ;
        Data{t+1,3} = CA;
        Data{t+1,4} = xA;
        save(filename, 'Data','-v7.3');
    end

    function merge_pops
        X = [xJ xA];
        if size(X,1) == size(unique(X, 'rows'),1)
        else
            [~, ind] = unique(X, 'rows', 'stable');
            X_uniq=X(ind,:);
            dup_ind = setdiff(1:size(X, 1), ind, 'stable')';
            X_dup=X(dup_ind,:);
            for i = 1:length(dup_ind)
                for j = 1:length(ind)
                    if X_dup(i,:) == X_uniq(j,:)
                        CA(ind(j)) = CA(ind(j)) + CA(dup_ind(i));
                        CJ(ind(j)) = CJ(ind(j)) + CJ(dup_ind(i));
                    else
                    end
                end
            end
            CA(dup_ind) = [];
            CJ(dup_ind) = [];
            X(dup_ind,:) = [];
            eJ1(dup_ind) = [];
            eJ2(dup_ind) = [];
            eA1(dup_ind) = [];
            eA2(dup_ind) = [];
        end
        xJ = X(:,1);
        xA = X(:,2);
    end
end