function cstr = setStateConstraint(value, knotpoint, index)
% Set State Constraint: Creates a structure describing a state constraint
%
%   Luke Drnach
%   February 28, 2020
cstr.value = value;
cstr.knotPoint = knotpoint;
cstr.stateIdx = index;
end