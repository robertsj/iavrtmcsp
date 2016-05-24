% SAMPLE INPUT FOR 1-D SLAB TRANSPORT w/ AUTOMATED VARIANCE REDUCTION
%  J. Roberts, 4/13/2010
clear
%profile on


% INPUT PARAMETER DEFINITIONS
%   The following provides descriptions of all possible user specified
%   variables for the code.  The sample problem specifications are defined
%   following the subsection of descriptions.
%
%   Problem Specification
%      xcm      = coarse mesh boundaries defining "cells"
%      numg     = number of energy groups
%      numm     = number of materials in cross-section array
%      xsec     = cross-section data for all materials
%                 has the format:
%                       mat 1  --> (T,g1) (A,g1) (S,g1->1) (S,g1->2) ...
%                                  (T,g2) (A,g2) (S,g2->1) (S,g2->2) ...
%                                  ...
%                       mat 2  --> ...
%                 (Isotropic scattering in the lab only!!!)
%      mt       = material assignment for each slab cell
%      src      = group wise source strength for each cell (note, this
%                 is renormalized in the code so that the integrated source
%                 distribution is unity)
%      det      = detector specification; this does two things: 1) it
%                 specifies which locations should get FOM's and 2) it
%                 defines the location (and allows energy-dependence)
%                 of the adjoint source if needed for VR. (The adjoint 
%                 source is the cross-section of the detector if all ones 
%                 are used in a specific region; otherwise, a differenct
%                 energy dependence is obtained (i.e. if only thermal
%                 reactions are of interest).
%      N        = number of histories

    xcm     = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20]/2;
    numg    = 1;
    numm    = 1;
    mt      = [1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1];    
    xsec    = [ 1.5      0.5      1.0 ];
               
    src     = [1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0];    
    det     = [0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1]; 
    
%   Variance Reduction Options
%      impcapt  = implicit capture [ 0=off, 1=on ]
%      srcbias  = source biasing [ 0=off, 1=on] (note, this requires that
%                 the user modify the appropriate portion of the source
%                 routine; instructions can be found there)
%      geosplt  = cell-based geometry splitting [ 0=off, 1=on ] (note, this
%                 uses the Sn-computed adjoint information)
%      stcadis  = weight windows and source biasing [ 0=off, 1=on ] (note,
%                 this uses the (standard) CADIS approach of Wagner, et al. 
%                 and relies on the adjoint information.  Weight windows 
%                 are specified on a grid independent of the actual 
%                 geometry. 
%      fwcadis  = weight windows and source biasing for global variance
%                 reduction [ 0=off, 1=on ] (note, this uses the
%                 forward-weighted CADIS approach.  Basically, the adjoint
%                 is weighted by a forward flux, and the end result can be
%                 better global statistics (i.e. for multiple or global
%                 tallies) then CADIS alone (which focuses more or less on
%                 a source-detector problem)
%      coopers  = weight windows and source biasing via a method similar
%                 to Cooper's method.  Essentially, the inverse foward flux
%                 is used in place of the adjoing.
%      wcut     = weight cutoff for implicit capture
%      wavg     = average weight for implicit capture (after roulette)
%      mcparts  = along with regular tallies, tally the UNWEIGHTED
%                 particles to estimate the population of monte carlo
%                 particles
%      CAUTION: Using more than one VR option might be bad bad news.


    N=1e5; 
    impcapt = 0; 
    srcbias = 0;
    geosplt = 0;
    stcadis = 1;
    fwcadis = 0;
    coopers = 0;

    wcut    = 1e-2;
    wavg    = 0.5;
    mcparts = 1; 

%  Discrete Ordinates Options (in addition to above)
%     xfm    = fine mesh interval for coarse divisions
%     ord    = number of ordinates (2,4,8, or 12)
%     maxit  = maximum iterations
%     maxerr = maximum relative pointwise error in phi

    xfm     = [ 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]*20;
    for o = 1:5
        if o == 1
            ord = 2;
        elseif o == 2
            ord = 4;
        elseif o == 3
            ord = 8;
        elseif o == 4
            ord = 16;
        else
            ord = 32;
        end
        
        maxit   = 100;
        maxerr  = 1e-6;
        
        input   =   struct( 'xcm',      xcm, ...
            'numslabs', length(xcm)-1, ...
            'numg',     numg, ...
            'numm',     numm, ...
            'xsec',     xsec, ...
            'mt',       mt, ...
            'src',      src, ...
            'det',      det, ...
            'N',        N, ...
            'impcapt',  impcapt, ...
            'srcbias',  srcbias, ...
            'geosplt',  geosplt, ...
            'stcadis',  stcadis, ...
            'fwcadis',  fwcadis, ...
            'coopers',  coopers, ...
            'wcut',     wcut, ...
            'wavg',     wavg, ...
            'mcparts',  mcparts, ...
            'xfm',      xfm, ...
            'ord',      ord, ...
            'maxit',    maxit, ...
            'maxerr',   maxerr);
        % call driver
        [out,snout] = mcdriver(input);
        a=1;
    end


