classdef ContactSolver
    
   methods (Abstract)
      [f, r] = solve(P, z);
      df = gradient(f, P, z, dP, dz);
   end
end