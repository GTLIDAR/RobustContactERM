function solverOptions = setSolverOptions(solverOptions, varargin)

% Create an input parser
parser = inputParser;
validPosScalar  = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validPosInteger = @(x) isnumeric(x) && isscalar(x) && (rem(x, 1) == 0) && (x > 0);
validScaleOpt = @(x) any(x == [0,1,2]);
% Add default values for the required fields
addParameter(parser, 'MajorFeasibilityTolerance', 1e-6, validPosScalar);
addParameter(parser, 'MajorOptimalityTolerance',1e-6, validPosScalar);
addParameter(parser, 'MinorFeasibilityTolerance',1e-6, validPosScalar);
addParameter(parser, 'ScaleOption',2, validScaleOpt);
addParameter(parser, 'MajorIterationsLimit',1000, validPosInteger);
addParameter(parser, 'IterationsLimit',10000, validPosInteger);
addParameter(parser, 'SuperbasicsLimit',500, validPosInteger);
addParameter(parser, 'ElasticWeight', 10^4, validPosScalar);
% Check the inputs
parse(parser, varargin{:});
% Return the defaults if necessary
if nargin == 0 || isempty(solverOptions)
    solverOptions = parser.Results;
else
    % Loop over the fields of the results and copy over only those which are
    % not defaults
    fields = fieldnames(parser.Results);
    for n = 1:length(fields)
       if ~any(strcmpi(fields{n}, parser.UsingDefaults))
          solverOptions.(fields{n}) = parser.Results.(fields{n}); 
       end
    end
end
end