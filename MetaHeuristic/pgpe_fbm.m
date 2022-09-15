decvar_init = [50; 0.723; 12];
fitness_init = -1.9813;
n_sample = 5;

sigma = [4, 1.5, 1.5]';

decvar_mid = decvar_init;
fitness_iter = fitness_init;
nEvo = 1;
while ~(nEvo >= 40 || min(fitness_iter) <= -1.995 || (nEvo >= 20 && min(fitness(end-19:end)) == max(fitness(end-19:end))))
    disp(['Number of Evolutions = ' num2str(nEvo)]);
    disp(['Fitness = ' num2str(fitness_iter(end))]);
    disp(['best decision Variables:' ]);
    disp(['N_freq = ' num2str(decvar_mid(1,end)) ]) ;
    disp(['Phasewarp = ' num2str(decvar_mid(2,end))]);
    disp(['T_wave = ' num2str(decvar_mid(3,end))]);
    
    d_dec1 = randn(3, n_sample).*sigma;
    d_dec2 = -d_dec1;
    decvars = decvar_mid + [d_dec1, d_dec2];
    fitness = match_random(decvars);
    d_dec_mid = sum((d_dec1 - d_dec2).*(fitness(1:n_sample)-fitness((n_sample+1):(2*n_sample))), 2)/sum(fitness);
    decvar_mid = decvar_mid - d_dec_mid;
    fitness_iter(end+1) = match_random(decvar_mid);
end

%%
function out = match_random(d)
    out = vecnorm(d, 1, 1);
end
