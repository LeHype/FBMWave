function fitness = match_random(decvars)
Bins = [-1E9:5E6:1E9]; % So far the best way i found was to check the input function
                       % and manually adjust the bins. 
use_DF = false;
use_derivative = true;
                       
                       
dt = 0.25;               % Time increment
T = 1000;                % Time Horizon
N_Seeds = 18;           % Number of different seeds used
persistent H1          % Histogramm of input Wave is the same everytime
persistent DH1         % Histogramm of Derivative
% persistent DF1         % Dominant Frequency
% iterate paralelly over all members
fitness = zeros(size(decvars,2),1);
fitness_abs = zeros(size(decvars,2),1);
fitness_der = zeros(size(decvars,2),1);
fitness_DF = zeros(size(decvars,2),1);
for i = 1: size(decvars,2)
    disp(i)
    time=[0:dt:T];
    N_freq = max(2,round(decvars(1,i),0));   %% Number of frequencies is integer and >1
    Phasewarp = round(abs(decvars(2,i)),0);   %% Phasewarp
    T_Wave = abs(decvars(3,i));      %% Period of dominant wave
%     AmpMax = abs(decvars(4,i));      %% Mean amplitude

    if isempty(H1)
        disp('beep boop')
    F1 = @(s,t) arrayfun(@(t) NewStochasticWave(t,'randomseed',s),t);
    seeds1 = [30:1:30+N_Seeds];
    S1 = zeros(length(seeds1),length(time));
    parfor j = 1:length(seeds1)
        S1(j,:) = F1(seeds1(j),time);
    end
%     S1.wait;
%     S1 = S1.fetchOutputs;
    S1 = reshape(S1',1,[]);
    h1 = histcounts(S1,Bins);
    H1 = h1;
    SS1 = S1;
    dS1 = (SS1(3:end)-SS1(1:end-2))/(2*dt);
    DH1 = histcounts(dS1,Bins);
%     DF1 = FindDominantFrequency(S1,1/dt);
    end    
    F2 =@(s,t) arrayfun(@(t) FBMStochasticWave(t,'Seed',s,'Phasewarp',Phasewarp,'N_freq',N_freq,'T_wave',T_Wave),t);  
    
    seeds2 = [1:1:1+N_Seeds];
    S2 = zeros(length(seeds2),length(time));
    parfor j = 1:length(seeds2)
        S2(j,:) = F2(seeds2(j),time);
    end
    
%     S2.wait;
%     S2out = S2.fetchOutputs;
    S2 = reshape(S2',1,[]);
    h2 = histcounts(S2,Bins);
    H2 = h2;
    SS2 = S2;
    dS2 = (SS2(3:end)-SS2(1:end-2))/(2*dt);
    DH2 = histcounts(dS2,Bins);
%     DF2 = FindDominantFrequency(S2,1/dt);
    
    fitness_der(i) = -sum((DH1- mean(DH1)).*(DH2- mean(DH2)))/sqrt(sum((DH1-mean(DH1)).^2).*sum((DH2-mean(DH2)).^2));
    
    fitness_abs(i) = -sum((H1- mean(H1)).*(H2- mean(H2)))/sqrt(sum((H1-mean(H1)).^2).*sum((H2-mean(H2)).^2));
    
%     fitness_DF(i)  = ((DF1-DF2)/DF1)^2;
    
end
fitness = fitness_abs;
if use_derivative

    fitness = fitness+fitness_der;
end 
% if use_DF
%     fitness = fitness +fitness_DF;
% end
end