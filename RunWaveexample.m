t = 100;    % How long to run the wave
Nfreq = 8; % How many frequencies to run



IndWaves = zeros(t,Nfreq); 
Waves = zeros(t,1);


for i = 1:t
    [wave, Amplitudes,Individual,Prediction] = FBMStochasticWave(i,'N_freq',Nfreq,'Seed',1,PhaseWarp=0.5,AmpWarp=0);
    IndWaves(i,:) = Individual;
    Waves(i) = wave;
    
end
plot(Waves)
%%
figure(1)
plot(Waves)
mean(abs(Waves))
%%
figure(2)
plot(IndWaves)
