function [total_cost,best_solution,total_violation ] = Start(sppnw, maximumIteration)

% Read sppnw file
[matrix_a, column_cost] = ReadInData(sppnw);

PopSize = 100;
iterator = 0;
isTerminated = false;
column_size = size(matrix_a,2);

% Call Initialize Function 
% Initialisation method (Algorithm 2, p341 in [1]) By using StochasticSetCover.m
[population, violation] = Initialize(matrix_a, column_cost, PopSize); 

% Call fitnessFunc Function 
fitness = fitnessFunc(population, matrix_a, column_cost); 

while isTerminated == false

    % Apply Stochastic ranking
    % StochasticRanking algorithm as described in [2] for constraint handling
    population = StochasticRanking(population, fitness, violation); 
    parent_val = PopSize * 0.3;
    parents = population(1:parent_val,:);
    offerspring = parents;

    % Apply Crossover
    index  = randi([1 parent_val],2, 1);
    cPoint =  randi([1 column_size]);
    offerspring(index(1,:), cPoint+1:end) = parents(index(2,:), cPoint+1:end);
    offerspring(index(2,:), cPoint+1:end) = parents(index(1,:), cPoint+1:end);
    
    % Apply Mutation
    for i=1:parent_val
        mutation = randi([1 column_size]);
        offerspring(i, mutation) = ~parents(i, mutation);
    end
    
    % Apply Heuristic improvement
    % Heuristic improvement Algorithm 1, p331 in [1]
    offerspring = HeuristicImprovement(offerspring, matrix_a, column_cost);
    
    population = [population; offerspring];   
    [fitness, total_cost, violation] = fitnessFunc(population, matrix_a, column_cost);
  
    [fitness, sorte_index] = sort(fitness);
    population = population(sorte_index,:);
    population = population(1:PopSize,:);
    
    iterator = iterator + 1;
    r = rem(iterator,100);
    
    if(r == 0)
        disp(['Please be patient, Iteration #: ' , num2str(iterator)]);
    end
    
    if(iterator >= maximumIteration)
        isTerminated = true;        
    end
end

best_solution = population(1,:);
[fitness, total_cost, total_violation]  = fitnessFunc(best_solution, matrix_a, column_cost);
disp('Process has finished. You can see the results:');
disp(char(10));
disp(['Minimum Cost: ', num2str(total_cost)]);
disp(['Violation: ', num2str(total_violation)]);
disp(['Solution: ', num2str(best_solution)]);

     
end

% Initialisation method (Algorithm 2, p341 in [1]) By using StochasticSetCover.m
function [population, V] = Initialize(matrix_coo, cCost, pSize)

D = size(matrix_coo, 2);
sz = pSize;
V = zeros(sz,1);
population = zeros(sz,D);
cost = zeros(sz,1);

for j = 1 : sz 
    [cost(j,1), population(j,:), V(j,1)] = StochasticSetCover(matrix_coo, cCost);
end

end

function [fitness, total_cost, tDegree]  = fitnessFunc(p, matrix_a, column_cost)
    pSize = size(p,1);
    total_cost = zeros(1, pSize);
    for k=1:pSize
       z = p(k,:);
       total_cost(1,k) = sum(z.*column_cost);
    end
    tDegree = zeros(1, pSize);
    for k=1:pSize
        z = p(k,:);
        h  =  (matrix_a*z'-1).^2;
        tDegree(1,k) = sum(h);
    end
    fitness = total_cost + 20000 * tDegree; 
end

% StochasticRanking algorithm as described in [2] for constraint handling
function [sort_population] = StochasticRanking(population, fit, p)
Pf = 0.45;
sort_population = population;
sort_fit = fit;
populationSize = size(population);

for m = 1 : populationSize
    isSwapped = false;
    for n = 1 : (populationSize - 1)
        u = rand(1);
        if (p(n) == p(n+1) == 0 || u<Pf)
            if (fit(n) > fit(n+1))
                SwapLists(sort_population,n,(n+1)); 
                SwapLists(sort_fit,n, (n+1));
                isSwapped = true;
            end
        else
            if p(n) > p(n+1)
                 SwapLists(sort_population,n, (n+1)); 
                 SwapLists(sort_fit, n, (n+1));
                isSwapped =true;
            end
        end        
        if isSwapped == false
            break;
        end
    end
end
end

function alteredList = SwapLists( List, i, j )
    alteredList = List(i);
    alteredList(i) = List(j);
    alteredList(j) = List(i);
end

% Heuristic improvement Algorithm 1, p331 in [1]
function freshPop = HeuristicImprovement(offspring, matrix_a, column_cost )

freshPop = offspring;
[C,Z] = size(freshPop);
nr = size(matrix_a,1);

for m = 1:C    
    ns = offspring(m,:);
    for n = 1:nr
        cr_m = find(matrix_a(n,:)==1);
        covered_ns = find(ns(1,:)==1);
        [q(n), Z] = size(intersect(cr_m, covered_ns));
    end 
    N_ns = ns;
    while sum(N_ns)>0
        covered_index = find(ns==1);
        rand = randi(length(covered_index));
        h = covered_index(rand);
        N_ns(h) = 0;
        for d = 1:nr
            if q(d) > 1
                ns(d) = 0;
                q(d) = w(h)-1;
            end
        end
    end    
    uncovered_index = find(q==0);
    [X,y] = size(uncovered_index);
    while y > 0
        rd = randi(length(uncovered_index));
        q(1,rd) = 0;
        h = uncovered_index(rd);
        cr_m = find(matrix_a(:,h)==1);        
        [minimum_cost, index] = min(column_cost(cr_m));
        index_arrayList = find(column_cost(cr_m) == minimum_cost);
        minimum_cost_tmp = length(index_arrayList);        
        if minimum_cost_tmp > 1
            index = index_arrayList(randi(minimum_cost_tmp));
        end
        h = cr_m(index);
        ns(1,h) = 1;  
    end    
    freshPop(m,:) = ns;
    
end

end

