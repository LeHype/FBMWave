t = 200; % How long to run the wave
Nfreq = 15; % How many frequencies to run



IndWaves = zeros(t,Nfreq); 
Waves = zeros(t,1);
for i = 1:t


[wave Individual] = SimplexStochasticWave(i,'N_freq',Nfreq,'Seed',2);
IndWaves(i,:) = Individual;
Waves(i) = wave;
end
% figure
% plot(IndWaves(:,2))
figure
plot(Waves)
%%
figure
for i = 1:Nfreq
    plot(IndWaves(:,i))
    hold on 
end