classdef ES_Individual < handle
    %ES_Individual This class defines an individual in an ES algorithm
    %
    %   This class is intended to create a data structure
    %   to represent an individual in an Evolution Strategy
    %   algorithm.
    %   Briefly, ES works by Recombination, Mutation and Selection of
    %   individuals in a set called Population in order to optimize a fitness
    %   function F.
    %
    %   Upon every iteration and until a certain stop condition is reached
    %   (e.g. number of iterations limit, mutation step size too low), ES
    %   will recombinate the population P to generate an Offspring. Then, P
    %   and its Offspring is mutated. Every individual is then evaluated
    %   with respect to F and the best are selected.
    %
    %   Every Individual is represented as a Cromossome, with 3 components:
    %   - x: a n-vector containing the variables to be evaluated
    %   - sigma: mutation step sizes, which may be a single value or a
    %   n-vector
    %   - alfa: rotation angles, a k-vector used for correlated mutation
    %
    %   where n is the dimension of the problem and k the possible number
    %   of pairs of variables (i.e. n*(n-1)/2)
    %
    %   In addition to that, to simplify further work, a property called
    %   f_eval was introduced. It corresponds to the fitness evaluation to
    %   a certain x
        
    properties
        x
        sigma = 1;
        alfa
        f_eval
    end % properties
    
    methods
        % Constructor
        function ind = ES_Individual(n,x_0,mutationStepSize,rotationAngles)
            if nargin == 0
                error('No input arguments provided')
            elseif nargin == 1 
                % If no boundaries or x_0 are provided, initialize them in the
                % interval [-100,100]
                ind.x = zeros(n,1);
                for i = 1:n
                    ind.x(i) = 100*(2*rand()-1);
                end % for
            end
            
            if nargin >= 2
                % Either a boundary or a x_0 was provided. First check
                % which one
                if isequal(size(x_0),[1 n]) % x_0
                    ind.x = x_0';
                elseif isequal(size(x_0),[n 1]) % x_0
                    ind.x = x_0;   
                elseif isequal(size(x_0),[1 1]) % boundary
                    ind.x = zeros(n,1);
                    for i = 1:n
                        ind.x(i) = x_0*(2*rand()-1);
                    end % for
                else 
                    error('Invalid input for x_0. It must be either a single value or a n-vector')
                end % if
            end % if
            
            if nargin >= 3 
                % mutationStepSize was provided
                ind.sigma = mutationStepSize;
            end
            
            if nargin >= 4
                % rotationAngles was provided
                ind.alfa = rotationAnlgles;
            end
        end % function Constructor
        
        % Mutation
        function mutate(ind,learningRate)
            % mutate sigma
            ind.sigma = ind.sigma * exp(learningRate*rand());
            % mutate x
            ind.x = ind.x + ind.sigma*rand();
        end % function mutate
        
        % Recombination
        function child = recombinate(p1,p2,method)
            n = max(size(p1.x));
            if method == 'avg'
                x_0 = (p1.x + p2.x)/2;
            elseif method == 'rnd'
                x_0 = zeros(n,1);
                for i = 1:n
                   if rand() >= 0.5
                       x_0(i) = p1.x(i);
                   else 
                       x_0(i) = p2.x(i);
                   end % if
                end % for
            end % if
            child = ES_Individual(n,x_0);
        end % function
        
    end % methods
    
end

