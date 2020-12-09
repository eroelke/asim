
function [aero, bijs] = get_n3_betas(aero, br21, br31, br41)
aero.m = [150 60 40 20];
aero.rcs(1) = 1;
aero.rcs(2) = sqrt( (aero.m(2) * aero.rcs(1)^2) / ...
    (aero.m(1) * br21) );
aero.rcs(3) = sqrt( (aero.m(3) * aero.rcs(1)^2) / ...
    (aero.m(1) * br31) );
aero.rcs(4) = sqrt( (aero.m(4) * aero.rcs(1)^2) / ...
    (aero.m(1) * br41) );

if (aero.rcs(4) > aero.rcs(3)) || (aero.rcs(3) > aero.rcs(2)) ... 
 || (aero.rcs(2) > aero.rcs(1))
    error('Impossible Geometry')
end
% 
aero.rn = 0.75/2;
aero.cds = [1.05 1.05 1.05 1.05];
aero.cls = [0 0 0 0];
% 
for x = 1:4
    b(x) = aero.m(x)/( aero.cds(1) * pi * aero.rcs(x)^2 );
end
br41 = b(4)/b(1);
br31 = b(3)/b(1);
br21 = b(2)/b(1);
% br32 = b(3)/b(2);
% br43 = b(4)/b(3);
% br42 = b(4)/b(2);

bijs = [num2str(br21) '_' num2str(br31) '_' num2str(br41)];

end
