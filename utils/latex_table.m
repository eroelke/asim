% function to output latex table code

function latex_table(outs, row_names, col_names)

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
col_headers = append(col_headers, ['\\textbf{' col_names{1} '}']);
for i = 2:n_cols
    col_headers = append(col_headers, [' & \\textbf{' col_names{i} '}']);
end
col_headers = append(col_headers, '\\\\ \\hline');

for i = 1:n_rows % row
    rows{i} = ['        ' row_names{i}];
    out = outs(i);
    
    col(1) = median(out.haf);
    col(2)= std(out.haf);
    col(3) = min(out.haf);
    col(4) = max(out.haf);
    %col(5) = length(find(out.haf < out.traj.alt(1,1)));
    %col(6) = length(find(out.haf < 0));
    col(5) = median(out.dv_circ);
    col(6) = std(out.dv_circ);

    for j = 1:(n_cols - 1) %cols
        rows{i} = append(rows{i}, [' & ' num2str(round(col(j),2))]);
    end
    rows{i} = append(rows{i}, '\\\\\n');
end



fprintf(['\n' ... 
'\\begin{table}[h!]\n' ...
'    \\centering\n' ...
'    \\begin{tabular}{|c|c|c|c|c|c|c|} \\hline\n' ...
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
names = {'Case','Median $h_a$','$1\\sigma h_a$','Min $h_a$', ... 
    'Max $h_a$','Median $\\Delta v_{circ}$', '$1\\sigma \\Delta v_{circ}$'};
end

function names = default_row_names()
names = {'DSF','DI','ECF','ECFH','ECF-Pert','ECFH-Pert'};
end

