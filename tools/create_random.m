% create random.txt

fid = fopen('random.txt','w+');

seeds = [2:1000];

for ii = 1:(length(seeds)-1)
    fprintf(fid,'%d\n',seeds(ii));
end
fprintf(fid,'%d',seeds(end));

fclose(fid);
        
