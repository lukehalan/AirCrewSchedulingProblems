function [total_cost, F, best_constraint] = StochasticSetCover(matrix_a, column_cost)

% number of rows
m = size(matrix_a,1);
% number of columns
n = size(matrix_a, 2);

% Initiate I to indicate all rows are not cover (I_i = 0 means the i_th row is not covered)
I = zeros(1,m);
% F is the solution, i.e., F_j=1 means the j_th column is selcted, F_j=0,
% otherwise
F = zeros(1,n);


% It terminates when all rows are covered, i.e., sum(I)=m
while sum(I)<m     
    % Find out which rows have not been covered
    uncovered_rows_idx = find(I==0);  
    % randomly select an uncovered row i
    i = uncovered_rows_idx(randi(length(uncovered_rows_idx)));
    % alpha_i is the indices of columns that cover row i
    alpha_i = find(matrix_a(i,:)==1);
    %  select column j \in \alpha_i which covers row i with minimum cost
    [mincost, idx] = min(column_cost(alpha_i));
    % However, there are multiple column with the same minimum cost
    % If we use min function in matlab, we will always selet the first column with the minimum cost
    idx_array = find(column_cost(alpha_i) == mincost);
    num_same_mincost = length(idx_array);
    % To prevent this problem, we randomly select one column if there are multiple column with the same minimum cost
    if num_same_mincost > 1
      idx = idx_array(randi(num_same_mincost));
    end
    j = alpha_i(idx);
    % Set column j as part of the solution
    F(1,j) = 1;
    % Set I to include the rows covered by column j
    I(1, matrix_a(:,j)==1) = 1;
end
total_cost = F*column_cost';
%disp(['The minimum cost found by the Stochastic Greedy algorithm is: ', num2str(total_cost)]);
%disp(['The solution found by the Stochastic Greedy algorithm is: ', num2str(F)]);
temp2 =  (matrix_a*F')';
best_constraint = sum(((temp2-1).^2)');
%disp(['The sum of volations of the constraints is: ', num2str(best_constraint)]);