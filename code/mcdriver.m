function [finout,snout] = mcdriver(in)
% function mcdriver(in)
%   This function takes as input the struct "in" containing all relevant
%   input parameters from the input file.  Driver then calls the
%   appropriate ''subroutines'' en route to the solution.

	disp('------------------------------------------------------')
    disp('--------- MULTI-GROUP SLAB MONTE CARLO WITH  ---------')
    disp('----- AUTOMATED ADJOINT-BASED VARIANCE REDUCTION -----')
    disp('------ Term Project for 22.106 by J. A. Roberts ------')
    disp('------------------------------------------------------')

% forward discrete ordinates -- used for comparison and possible vr
snout = sndriver(in);
disp('...finished forward SN')

% unbiased source distribution -- takes geometry and makes source densities
srcout = srcdriver(in);
disp('...finished source driver')

% variance reduction -- computes any required vr parameters
[vrout,snout] = vrparm(in,snout,srcout);
disp('...finished vr parameter production')

% play the game -- perform the monte carlo simulation
out = mcslab1d(in,srcout,vrout);
disp('...finished the game and now printing...')

% plots and such
finout = out1d(in,out,snout,srcout,vrout);


end