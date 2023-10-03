xAnalisis = [727.778,883.333,922.222,883.333,1000,1000,922.222];
yAnalisis = [1000,850,1000,900,850,750,800];

xLen = length(xAnalisis);
yLen = length(yAnalisis);

for i = 1:xLen
    load("./../Fase2/nModes/Waveguide"+string(xAnalisis(i)*1000)+"_"+string(yAnalisis(i)*1000)+"_1580.mat");
    disp(nmodes);
end