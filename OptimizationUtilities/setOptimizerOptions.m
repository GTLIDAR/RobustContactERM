function optimOptions = setOptimizerOptions(optimOptions, varargin)

% Luke Drnach
% February 28, 2020

% Check for an options subfield
if nargin == 0 || isempty(optimOptions) || ~isfield(optimOptions, 'options')
    [options, unmatched] = checkOptimizerOptions([], varargin{:});
else
    [options, unmatched] = checkOptimizerOptions(optimOptions.options, varargin{:});   
end
% Check the remaining input arguments
% Create an input parser
parser = inputParser;
validPosInteger = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (rem(x, 1) == 0);
validDurations = @(x) isnumeric(x) && all(x>0) && length(x) <= 2;

% Add default values for the parameter fields
addParameter(parser, 'nPoints', 101, validPosInteger);
addParameter(parser, 'duration',[1,1], validDurations);
addParameter(parser, 'display',false, @(x)islogical(x));


% Parse the unmatched inputs
parse(parser, unmatched{:});
if nargin == 0 || isempty(optimOptions)
    optimOptions = parser.Results;
else
    fields = fieldnames(parser.Results);
    for n = 1:length(fields)
        if ~any(strcmpi(fields{n}, parser.UsingDefaults))
            optimOptions.(fields{n}) = parser.Results.(fields{n});
        end
    end
end
% Add the additional options to the output structure
optimOptions.options = options;
end

function [options, unmatched] = checkOptimizerOptions(options, varargin)
% Helper function

% Luke Drnach
% February 28, 2020

parser = inputParser;
parser.KeepUnmatched = true;
% Add the parameter fields
if nargin == 0 || isempty(options)
    options = [];
end

isNonnegScalar = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
isPosScalar = @(x) isnumeric(x) && isscalar(x) && (x > 0);
isPosInteger = @(x) isPosScalar(x) && (rem(x, 1) == 0);
inrange= @(x, lb, ub) (isPosInteger(x) || x == 0) && (x >= lb) && (x <= ub);

% Set default values for contact-implicit trajectory optimization
addParameter(parser, 'integration_method', RobustContactImplicitTrajectoryOptimizer.BACKWARD_EULER, @(x) inrange(x, 1, 4));
addParameter(parser, 'uncertainty_source', RobustContactImplicitTrajectoryOptimizer.NO_UNCERTAINTY, @(x) inrange(x, 0, 3));
addParameter(parser, 'time_option', 1, @(x) inrange(x, 1, 2));
addParameter(parser, 'time_constraints', RobustContactImplicitTrajectoryOptimizer.BOUNDTIME, @(x) inrange(x, 1, 3));
addParameter(parser, 'nlcc_mode',2, @(x) inrange(x, 1, 5));
addParameter(parser, 'relax_cost',1, isPosScalar);
addParameter(parser,'compl_slack',0, isNonnegScalar);
% Set the default values for Robust-Specific options
addParameter(parser, 'distribution', RobustContactImplicitTrajectoryOptimizer.GAUSSIAN, @(x) inrange(x, 1, 2));
addParameter(parser, 'heightVariance', 1, isPosScalar);
addParameter(parser, 'frictionVariance', 1, isPosScalar);
addParameter(parser, 'ermScaling', RobustContactImplicitTrajectoryOptimizer.NOSCALE, @(x) inrange(x, 1, 3));
addParameter(parser, 'contactCostMultiplier',0 , isNonnegScalar);
addParameter(parser, 'ermMode',1,@(x) inrange(x, 1, 4));
addParameter(parser, 'distanceScaling',1, isPosScalar);
addParameter(parser, 'ermFrictionBias',0, isNonnegScalar);
addParameter(parser, 'frictionCostMultiplier', 1, isPosScalar);
addParameter(parser, 'distanceCostMultiplier', 1, isPosScalar);
% Parse the inputs
parse(parser, varargin{:});
% Get the results from the input parser
if nargin == 0 || isempty(options)
   options = parser.Results;
else
    fields = fieldnames(parser.Results);
    for n = 1:length(fields)
       if ~any(strcmpi(fields{n}, parser.UsingDefaults))
          options.(fields{n}) = parser.Results.(fields{n}); 
       end
    end
end
% Also return unmatched results
unmatchedFields = fieldnames(parser.Unmatched);
unmatched = cell(1, 2*length(unmatchedFields));
for n = 1:length(unmatchedFields)
    unmatched{2*n-1} = unmatchedFields{n};
    unmatched{2*n} = parser.Unmatched.(unmatchedFields{n});
end
end