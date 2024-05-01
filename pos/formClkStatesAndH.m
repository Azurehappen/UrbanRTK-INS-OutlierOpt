function [H_clk,x_clk] = formClkStatesAndH(num_sv)
% Create clock offset matrix for multi-GNSS
total = sum(num_sv);
num_sys = sum(num_sv~=0);
H_clk = zeros(total,num_sys);
row_ind = 0;
col_ind = 1;
for i = 1:length(num_sv)
    if num_sv(i) == 0
        continue;
    end
    H_clk(row_ind+1:row_ind+num_sv(i),col_ind)=1;
    row_ind = row_ind + num_sv(i);
    col_ind = col_ind + 1;
end
x_clk = zeros(num_sys,1);

end