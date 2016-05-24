function [vrout,snout] = vrparm(in,snout,srcout)
% function vrout = vrparm(in,snout)
%   This function generates any needed variance reduction parameters
%   required for the adjoint-based approaches (i.e. geometry splitting, or
%   either of the CADIS approaches.  Only one can be used!
tic

vrout.imp = 0;
% For reference, review the basic ideas of each variance reduction method
% option that uses the adjoint.
%
% Geometry Splitting 
%     Reference:  F.B. Brown, "Fundamentals of Monte Carlo Particle
%                 Transport", LA-UR-05-4983
%     In geometry splitting, an importance is assigned to each cell.  When
%     a neutron leaves cell a with importance ia, and goes to cell b with
%     importance ib, either ia<ib or ia>ib (or possibly, though unlikely,
%     ia=ib).  Then we let
%                  r = ib/ia   
%                  if r > 1
%                    n = floor(r)
%                    if n > 1
%                       split into n particles of w = w/n
%                       and each particle has the same (x,mu,g,...)
%                    end
%                  else % r < 1
%                    if rand < r
%                      w = w/r
%                    else
%                      kill particle
%                    end
%                  end
%     Note, it may or may not be wise to implement a weight cutoff with
%     splitting.
%
%  Parameters:
%     imp   = cell energy-dependent importances (just average adjoint
%             scalar flux over cell)
if in.geosplt==1

    % compute the adjoint
    [snout.phiA,snout.psiA,snout.phiAplot,snout.phiAavg,snout.xww] = ...
            sn_one_d(in,1);

    % the importance is just the average adjoint flux i.e. cell importance
    vrout.imp = snout.phiAavg;

% CADIS (Consistent Adjoint-Driven Importance Sampling
%     Reference: J.C. Wagner and A. Haghighat, "Automated Variance
%                Reduction of Monte Carlo Shielding Calculations Using the
%                Discrete Ordinates Adjoint Function", Nuclear Science and
%                Engineering, 128, 186-208 (1998)
%     The CADIS method does two things: it biases the transport of
%     particles via weight windows, and it biases the source so that
%     particles are born with weights withing their local weight window,
%     i.e. with consistent weights.  This consistency makes the method much
%     more efficient.  Note, while the underlying adjoint can theoretically
%     be a global one, CADIS works best for simple source-detector
%     problems.
%
%     Source Biasing
%          The source energy and position are sample from a biased source
%          distribution defined as:
%               q*(r,E)  =  [ phiA(r,E)*q(r,E) ] / < phiA(r,E)*q(r,E) >
%                        =  [ phiA(r,E)*q(r,E) ] / R
%          where "< >" means integration over all space and energies.  The
%          denominator is the total detector response.  Because the source
%          particles are biased, their weights are defined to be:
%                W(r,E)  =  R/phiA(r,E)
%          Note, "phiA" is the adjoint scalar flux.
%     Transport Biasing
%          First, recall the idea of weight windows.  For a particular cell
%          i, upper and low weight windows (wU and wL) are defined.  If a
%          particle of weight w enters cell i, we have three possibilities:
%          1) w < wL, 2) w>wU, or 3) wL<w<wU.  In the third case, nothing
%          happens, but case 1) leads to roulette and 2) leads to
%          splitting.  Here's the algorithm:
%               if w > wU
%                 n = min(mxspln, 1+wgt/wU) <- mxspln is max split ratio
%                 w = w/n
%                 bank n-1 copies, and follow the current particle
%               elseif w < wL
%                 P = max(1/mxspln, w/wavg) <= wavg is avg w 
%                 if rand < P
%                   w = w/P
%                 else
%                   kill particle
%                 end
%               end
%          Now to choose the wU and wL values we use CADIS, and so we have:
%               wL(r,E) =  [R/phiA(r,E)] * [2/(cU + 1)] 
%          where cU defines the width such that cU=wU/wL.  The default cU
%          in MCNP is 5.  We use that here, too.
elseif (in.stcadis+in.fwcadis+in.coopers) > 0
    
    % SOURCE BIASING
    
    %  step 1.  reconstruct unbiased source density mapped onto fine grid
    %     Recall that P is the cell/energy discrete source density
    %     in the form [ p(c1,g1)  p(c2,g1) ...
    %                   p(c1,g2) ...
    %     with corresponding cell bounds given by sx(ci,1) and sx(ci,2).
    %     What we need is a new discrete density corresponding to the same
    %     source mapped onto the fine grid defined by in.xfm.  Within a
    %     (coarse) cell, the source is uniformly distributed, and so should
    %     it be in a fine cell; hence, the probabilities for the fine cells
    %     will be just the original (coarse) probability divided by the
    %     number of fine cells, and the total number of columns of the new
    %     (still unbiased) source vector will be the total number of fine
    %     meshes in all source cells.  Hence,
    Pf  = zeros(in.numg, sum(in.xfm(srcout.srcloc)) );
    sxf = zeros(sum(in.xfm(srcout.srcloc)), 2);
    j = 0;
    for i = 1:length(srcout.srcloc)
        ii = srcout.srcloc(i);
        for g = 1:in.numg
            Pf( g, (j+1):(j+in.xfm(ii)) )  = ...
                srcout.P(g,i)/in.xfm(ii);
        end
        dx = ( srcout.sx(i,2)-srcout.sx(i,1) ) / in.xfm(ii);
        for k = 1:in.xfm
            sxf(j+k,1) = (k-1)*dx + srcout.sx(i,1);
            sxf(j+k,2) = dx+sxf(j+k,1);
        end
        j = j+in.xfm(ii);
    end
    
    %  step 2.  bias the resulting density and produce weights

    % compute the adjoint
    if in.stcadis == 1
        [snout.phiA,snout.psiA,snout.phiAplot,snout.phiAavg] =...
            sn_one_d(in,1);
    elseif in.fwcadis == 1
        [snout vrout] = fwcadis(in,snout,srcout,vrout);
    elseif in.coopers == 1
        %[snout.phiA,snout.psiA,snout.phiAplot,snout.phiAavg] =...
        %    sn_one_d(in,1);        
        snout.phiA = 1./snout.phiF;
    end
    
    pfbias = zeros(size(Pf)); % temporary
    weight = pfbias;
    for i = 1:srcout.numsrc
        % pick out the adjoint corresponding to this source location
        if srcout.srcloc(i) > 1
            step = sum( in.xfm(1:srcout.srcloc(i)-1) );
        else
            step = 0; step2 = 0;
        end
        pfbias(:,step2+1:step2+in.xfm(srcout.srcloc(i))) = ...
            ( snout.phiA(step+1:step+in.xfm(i),:).*...
              Pf(:,step2+1:step2+in.xfm(srcout.srcloc(i)))' )';
        weight(:,step2+1:step2+in.xfm(srcout.srcloc(i))) = ...
                1 ./ snout.phiA(step+1:step+in.xfm(i),:)';
        step2 = step2 + in.xfm(srcout.srcloc(i));
    end
    R      = sum(sum(pfbias)); % R is just the sum of all the q*phi's
    Pf     = pfbias / R;       % and we renormalize our source density
    Pfx    = sum(Pf);
    weight = weight * R;       % and our weights.

    vrout.Pf     = Pf;
    if (in.numg>1)
        vrout.Pxf = sum(Pf);
    else
        vrout.Pxf = Pf;
    end
    vrout.sxf    = sxf;
    vrout.weight = weight;
    vrout.numsrc = length(vrout.Pxf);
    
    % TRANSPORT BIASING (...remarkably more straightforward!)
    %   wL(r,E) =  [R/phiA(r,E)] * [2/(cU + 1)]
    vrout.xww(:,1) = snout.xww(1:end-1)';
    vrout.xww(:,2) = snout.xww(2:end)';
    vrout.cU = 5.0; % default ratio of upper and lower bounds
    vrout.wL = (2*R/(vrout.cU+1))./snout.phiA;
    
end

vrout.t = toc; % total vr parameter generation time
end

function [snout vrout] = fwcadis(in,snout,srcout,vrout)
% FW-CADIS (Forward-Weighted CADIS)
%     Reference: J.C. Wagner et al. "Forward-Weighted Cadis Method for
%                Variance Reduction of Monte Carlo Calculations of 
%                Distributions and Multiple Localized Quantities",
%                International Conference on Mathematics, Computational
%                Methods, and Reactor Physics (M&C 2009)
%          The basic idea of CADIS remains, but a minor modification is
%          made.  Instead of the original adjoint source, a forward
%          flux-weighted adjoint source is used.  For example, if our
%          desired "response" is uniform statistical error throughout the
%          problem, then (see the reference) it can be shown that we want
%                             Qadj(r,E) = 1 / phiF(r,E)
%          Of course, phiF is exactly what we aim to calculate!  But here, 
%          of course, we use a crude SN approximation for it.
%
%          Because the only difference between CADIS and FW-CADIS is the 
%          adjoint source used to generate the adjoint fluxes, FW-CADIS is
%          implemented here essentially as an optional subroutine called
%          from within the CADIS routine above.
    
    % first build the coarse adjoint distribution
    if in.numg>1
        srcloc = find(sum(in.det)>0); % picks out cell numbers w/ src
    else
        srcloc = find(in.det>0);
    end
    denom = 0;
    for i = 1:length(srcloc)
        V(i) = in.xcm( srcloc(i)+1)-in.xcm( srcloc(i) );
        sx(i,1)=in.xcm( srcloc(i) );  % lower bound
        sx(i,2)=in.xcm( srcloc(i)+1); % upper bound
        denom = denom + V(i)*sum( in.det(:,srcloc(i)));  
    end
    for i = 1:length(srcloc)
        P(1:in.numg,i)=V(i)*in.det(1:in.numg, srcloc(i))/denom;    
    end
    if in.numg > 1
        Px=sum(P);
    else
        Px = P;
    end
    numsrc=length(srcloc);
    % then build the fine adjoint distribution
    Pf  = zeros(in.numg, sum(in.xfm(srcloc)) );
    sxf = zeros(sum(in.xfm(srcloc)), 2);
    j = 0;
    for i = 1:length(srcloc)
        ii = srcloc(i);
        for g = 1:in.numg
            Pf( g, (j+1):(j+in.xfm(ii)) ) = P(g,i)/in.xfm(ii);
        end
        dx = ( sx(i,2)-sx(i,1) ) / in.xfm(ii);
        for k = 1:in.xfm
            sxf(j+k,1) = (k-1)*dx + sx(i,1);
            sxf(j+k,2) = dx + sxf(j+k,1);
        end
        j = j+in.xfm(ii);
    end
    % now weight the adjoint distribution with the forward flux
    qad = zeros(size(snout.phiF));
    step = 0; step2 = 0;
    for i = 1:numsrc
        % pick out the adjoint corresponding to this source location
        if srcloc(i) > 1
            step = sum( in.xfm(1:srcloc(i)-1) );
        else
            step = 0; step2 = 0;
        end
        qad(step+1:step+in.xfm(i),:) = ...
          ( Pf(:,step2+1:step2+in.xfm(srcloc(i)))' ./ ...
            snout.phiF(step+1:step+in.xfm(i),:) );
        step2 = step2 + in.xfm(srcloc(i));
    end
    % now calculate the adjoint using the weighted adjoint source
    [snout.phiA,snout.psiA,snout.phiAplot,snout.phiAavg,snout.xww] = ...
            sn_one_d(in,1,qad);

end
