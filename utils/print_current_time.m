% print_current_time.m
%   prints the current (12hr) time
% 
% Written by: Evan Roelke
% 
function print_current_time()

t = datetime;
if (t.Hour > 12)
    fprintf('%2.0f:%2.0f:%2.0f PM.\n',t.Hour - 12, t.Minute, t.Second);
else
    fprintf('%2.0f:%2.0f:%2.0f AM.\n',t.Hour, t.Minute, t.Second);
end
end