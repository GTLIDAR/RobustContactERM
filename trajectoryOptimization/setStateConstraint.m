function cstr = setStateConstraint(value, knotpoint, index)
% Set State Constraint: Creates a structure describing a state constraint.
%   This function is used in conjunction with setOptimizerOptions to store
%   the optimization options and state constraints between optimization
%   runs.
%
%   Arguments:
%       value: a vector of values for the state constraint
%       knotpoint: the index of the knot point at which to apply the
%       constraint
%       index: the state variable indices subjected to the constraint
%
%   Luke Drnach
%   February 28, 2020
cstr.value = value;
cstr.knotPoint = knotpoint;
cstr.stateIdx = index;
end