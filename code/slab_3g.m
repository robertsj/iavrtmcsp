% 3-GROUP SLAB PROBLEM
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
    numg    = 3;
    numm    = 1;
    mt      = [1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1];
    
%    if(numg==3)
    xsec    =  [   1.0  0.3  0.3  0.2  0.2
                   1.5  0.9  0.0  0.4  0.2
                   2.0  0.5  0.0  0.0  1.5  ];
    
    src     = [  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  
                 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ];    
                 
    det     = [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 ];
%     else       
%     xsec    = [ 1.5      0.5      1.0 ];
%                
%     src     = [1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0];    

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

for m = 1:10

    if m == 1 % analog
        N=4e7;impcapt = 0; srcbias = 0;geosplt = 0;stcadis = 0;fwcadis = 0;coopers = 0;
        det = [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 ];
    elseif m==2 % implicit capture
        N=1e7;impcapt = 1; srcbias = 0;geosplt = 0;stcadis = 0;fwcadis = 0;coopers = 0;
        det = [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 ];
    elseif m==3 % geometry splitting
        N=1e6;impcapt = 0; srcbias = 0;geosplt = 0;stcadis = 1;fwcadis = 0;coopers = 0;
        det = [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 ];
    elseif m==4 % st cadis, slab 7
        N=1e5;impcapt = 0; srcbias = 0;geosplt = 0;stcadis = 1;fwcadis = 0;coopers = 0;
        det = [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  1 ];
    elseif m==5 % st cadis, slab 14
        N=1e5;impcapt = 0; srcbias = 0;geosplt = 0;stcadis = 1;fwcadis = 0;coopers = 0;
        det = [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0 ];
    elseif m==6 % st cadis, slab 20
        N=1e6;impcapt = 0; srcbias = 0;geosplt = 0;stcadis = 1;fwcadis = 0;coopers = 0;
        det = [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 ];
    elseif m==7 % st cadis, three
        N=1e6;impcapt = 0; srcbias = 0;geosplt = 0;stcadis = 1;fwcadis = 0;coopers = 0;
        det = [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  0  1  0  0  0  0  0  0  1  0  0  0  0  0  1 ];
    elseif m==8 % fw cadis, three
        N=1e6;impcapt = 0; srcbias = 0;geosplt = 0;stcadis = 0;fwcadis = 1;coopers = 0;
        det = [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  0  1  0  0  0  0  0  0  1  0  0  0  0  0  1 ];
    elseif m==9 % fw cadis, all
        N=1e6;impcapt = 0; srcbias = 0;geosplt = 0;stcadis = 0;fwcadis = 1;coopers = 0;
        det = [0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
               0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
               1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1];
    elseif m==10 % coopers
        N=1e6;impcapt = 0; srcbias = 0;geosplt = 0;stcadis = 0;fwcadis = 0;coopers = 1;
        det = [1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1
               1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1
               1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1];
    end
    wcut    = 1e-2;
    wavg    = 0.5;
    mcparts = 1; 

%  Discrete Ordinates Options (in addition to above)
%     xfm    = fine mesh interval for coarse divisions
%     ord    = number of ordinates (2,4,8, or 12)
%     maxit  = maximum iterations
%     maxerr = maximum relative pointwise error in phi

    xfm     = [ 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]*10;
    ord     = 32;
    maxit   = 100;
    maxerr  = 1e-6;    
     
% DO NOT MODIFY BELOW
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

if m == 1
    save 'slab_3g_analog'
elseif m == 2
    save 'slab_3g_impcapt'
elseif m == 3
    save 'slab_3g_geomsplt'
elseif m == 4
    save 'slab_3g_stcadis_cell_7'
elseif m == 5
    save 'slab_3g_stcadis_cell_14'
elseif m == 6
    save 'slab_3g_stcadis_cell_20'
elseif m == 7
    save 'slab_3g_stcadis_all_three'
elseif m == 8
    save 'slab_3g_fwcadis_all_three' 
elseif m == 9
    save 'slab_3g_fwcadis_full'
elseif m == 10
    save 'slab_3g_coopers'
end

end
%profile viewer


