


function A = GetPlanetaryWaveAmplitude(z,param,lat)

params = {'density','pressure','temperature','windE','windN','windU'};


if (z >= 45 && z <= 65)
    switch param
        case 1 %density
            
        case 2 %pressure
        case 3 %temp
        case 4 %windE
        case 5 %windN
        case 6 %windU
            
        otherwise
            A = 0;
    end
end




%{

Planetary Scale Waves
        low_lat  mid_lat  high_lat
Au      10       9.7      8.1
Av      1.1      6.3      9.8
Ap      0.4      3.6      6.3
Arho    0.6      2.7      6.9
AT      0.4      0.9      1.7

%}


end