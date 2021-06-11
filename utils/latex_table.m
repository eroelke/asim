% function to output latex table code

function latex_table(outs, row_names, col_names)

if (isstruct(outs))
    if isfield(outs, 'out_hyb')
        outs = [outs.out_k, outs.out_di, outs.out_ecf, outs.out_hyb, outs.out_ecf_p, outs.out_hyb_p];
    else
        outs = [outs.out_k, outs.out_di, outs.out_ecf, outs.out_ecfh, outs.out_ecf_p, outs.out_ecfh_p];
    end
end

switch (nargin)
%     case 0 % all defaults
%         if (exist('out_dsf', 'var') == 0)
%             outs = [out_k, out_di, out_ecf, out_hyb, out_ecf_p, out_hyb_p];
%         else
%             outs = [out_dsf, out_di, out_ecf, out_hyb, out_ecf_p, out_hyb_p];
%         end
%         row_names = {'DSF','DI','ECF','ECFH','ECF-Pert','ECFH-Pert'};
%         col_names = {'Case','Mean $h_a$','$1\\sigma h_a$','Min $h_a$','Max $h_a$'};
    case 1 % default row, col labels
        row_names = default_row_names();
        col_names = default_col_names();
    case 2 % default col names
        col_names = default_col_names();
    case 3
        % default nothing
    otherwise
        % default nothing
end

n_rows = length(row_names);
n_cols = length(col_names);

% rows = ['        DSF & x & x & x & x \\\\\n' ...
% '        DI & x & x & x & x \\\\\n' ...
% '        ECF & x & x & x & x \\\\\n' ...
% '        ECFH & x & x & x & x \\\\\n' ...
% '        ECF-Pert & x & x & x & x \\\\\n' ...
% '        ECFH-Pert & x & x & x & x \\\\' ];

col_headers = '        ';
col_headers = append(col_headers, [col_names{1}]);
for i = 2:n_cols
    col_headers = append(col_headers, [' & ' col_names{i}]);
end
col_headers = append(col_headers, '\\\\ \\hline');

for i = 1:n_rows % row
    rows{i} = ['        ' row_names{i}];
    out = outs(i);
    
%     for j = 1:n_cols - 1
%        switch col_headers{j}
%            case  
%            
%        end
%     end
    
    nCrashes = length(find(out.haf < out.traj.alt(1,1)));
    nEscapes = length(find(out.haf < 0));
    
    col(1) = mean(out.haf_err);
    col(2) = std(out.haf_err);
%     col(3) = col(1) - 3*col(2);
%     col(4) = col(1) + 3*col(2);
    col(3) = iqr(out.haf_err);
    col(4) = prctile(out.haf_err, 95);
    
    col(5) = min(out.haf);
    if (col(5) <= 51)
        col(5) = 0;
    end
    col(6) = max(out.haf);
    col(7) = mean(out.dv_circ);
    col(8) = std(out.dv_circ);
    col(9) = nanmean(out.tjr(:,1));
    col(10) = nanstd(out.tjr(:,1));

    for j = 1:(n_cols - 1) %cols
        if (j > 6)
            rows{i} = append(rows{i}, [' & ' num2str(round(col(j),3))]);
        else
            rows{i} = append(rows{i}, [' & ' num2str(round(col(j),1))]);
        end
    end
    rows{i} = append(rows{i}, '\\\\\n');
end



fprintf(['\n' ... 
'\\begin{table}[h!]\n' ...
'    \\centering\n' ...
'    \\begin{tabular}{|' repmat('c|',[1, n_cols])  '} \\hline\n' ...
col_headers '\n']);

for i = 1:n_rows
    fprintf(rows{i});
end
    
fprintf([ ... 
'        \\hline\n' ...
'    \\end{tabular}\n' ...
'    \\caption{caption}\n' ...
'    \\label{tab:label}\n' ...
'\\end{table}\n' ...
]);



end

function names = default_col_names()
% 
% names = {'\\textbf{Case}' , ...                 % 1
%     '$\\Delta h_a \\: \\mu$' , ...              % 2
%     '$\\sigma \\: \\Delta h_a$' , ...           % 3
%     '$\\Delta h_a \\: \\mu-3\\sigma$' , ...     % 4
%     '$\\Delta h_a \\: \\mu+3\\sigma$' , ...     % 5
%     'Min $h_a$' , ...                           % 6
%     'Max $h_a$' , ...                           % 7
%     '$\\Delta v_{c} \\: \\mu$' , ...            % 8
%     '$\\sigma \\: \\Delta v_{c}$' , ...         % 9
%     '$t_j/t_f \\: \\mu$' , ...                  % 10
%     '$\\sigma \\: t_j/t_f$' ...                 % 11
%     };

names = {'\\textbf{Case}' , ...                 % 1
    '$\\Delta h_a \\: \\mu$' , ...              % 2
    '$\\sigma \\: \\Delta h_a$' , ...           % 3
    '$\\Delta h_a$ IQR' , ...     % 4
    '\\begin{tabular}{@{}c@{}}$\\Delta h_a$ \\\\ 95th Percentile\\end{tabular}' , ...     % 5
    'Min $h_a$' , ...                           % 6
    'Max $h_a$' , ...                           % 7
    '$\\Delta v_{c} \\: \\mu$' , ...            % 8
    '$\\sigma \\: \\Delta v_{c}$' , ...         % 9
    '$t_j/t_f \\: \\mu$' , ...                  % 10
    '$\\sigma \\: t_j/t_f$' ...                 % 11
    };
end

function names = default_row_names()
names = {'DSF','DI','ECF','ECFH','ECF-P','ECFH-P'};
end

