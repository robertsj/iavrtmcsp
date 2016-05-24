function snout = sndriver(in)
% function snout = sndriver(in)
%    This function takes as input the struct "in" which is defined in the
%    sample input.  The function calls the local function sn_one_d.
%    J. Roberts, 4/13/2010

% compute the forward quantities
tic

[snout.phiF,snout.psiF,snout.phiFplot,snout.phiFavg,snout.xww,snout.x] = ...
    sn_one_d(in,0);

snout.snt = toc;

end

