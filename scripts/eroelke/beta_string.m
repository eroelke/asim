function s = beta_string(n, index)

prefix = '\beta_{i+1}/\beta_i = \{';
switch (n)
    case 1
        switch (index)
            case 1
                s = '\beta_2/\beta_1 = 10';
            otherwise
                s = '';
        end
    case 2
        switch (index)
            case 1
                s = [prefix '3,10\}'];
            case 2
                s = [prefix '5,10\}'];
            case 3
                s = [prefix '7,10\}'];
            case 4
                s = [prefix '9,10\}'];
            otherwise
                s = '';
        end
    case 3
        switch (index)
            case 1
                s = [prefix '3,5,10\}'];
            case 2
                s = [prefix '5,7,10\}'];
            case 3
                s = [prefix '7,9,10\}'];
            case 4
                s = [prefix '9,9.5,10\}'];
            otherwise
                s = '';
        end
    case 4
        % to do
    otherwise
        s = '';
end

end