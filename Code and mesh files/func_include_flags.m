%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======================= DEFINE GLOBAL VARIABLES =========================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FEM MODEL PARAMETERS 

% NONLOCAL GRADIENT PARAMETER 
% g = lc^2/2, lc = characteristic length
global g

% ADAPTIVE LOAD AND PLOTTING PARAMETERS
global min_iter max_iter max_accept_iter dlfactor_incr_threshold increment_plot_threshold loadfactor_plot_threshold

% SOLVER SCHEME
% SolverID: 1 - Local, 2 - Nonlocal Gradient, 3 - Nonlocal Integral
% TangentID: 1 - Analytical, 2 - Numerical
global SolverID TangentID

global newDecomp L U P CA_inverse CA_inverseB lenA lenD
