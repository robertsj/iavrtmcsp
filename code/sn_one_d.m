function [phi,psi,phiPLOT,phiAVG,xe,xa] = sn_one_d(in,adj,qad)

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% This function solves the 1-D multigroup SN equations given input.  It's !
% output is the forward (or adjoint) angular fluxes at the cell edges and !
% the cell-centered scalar flux.  It has been verified against PARTISN for!
% some simple problems.
% ** last modified by J. Roberts, 4/11/2010
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% 1) struct "in" contains (among other things)
%     numg   = number of energy groups
%     numm   = number of materials
%     xcm    = coarse divisions
%     xfm    = fine mesh interval for coarse divisions
%     mt     = material assignment for each coarse division
%     xsec   = cross-sections in the form
%              (mat1/g1) sigTOT   sigA  sigSg1->g1  sigSg1->g2 ...
%                 ...g2) sigTOT   sigA  sigSg2->g1  ....
%              (mat2/g1) ...
%     src    = volumetric (isotropic) source by coarse mesh and energy
%              (for adj, src is the cross-section for the source region)
%     det    = spatial location of detector
%     ord    = number of ordinates (2,4,8, or 12)
%     maxit  = maximum iterations
%     maxerr = maximum relative pointwise error in phi
% 2)  adj    = 0-->forward, 1-->adjoint
% 3)  qad    = if present, it is a fine mesh adjoint source
% working variables
%     dx     = fine mesh divisions
%     mtt    = material assignment for fine mesh cells
%     n      = number of fine mesh cells

numg = in.numg;
ord  = in.ord;
Q = zeros(  sum(in.xfm), ord, numg );

bound       =   [ord/2+1 ord; 1 ord/2];
gbound      =   [ 1  numg];
git         =   1;
S           =   in.src;
if adj == 1
    git         = -1;
    gbound      = [numg 1];
    bound       = flipud(bound);
    % define S to be the group-wise TOTAL cross section in the detector
    S = in.det;
%     for i = 1:length(in.xfm) % number of course meshes
%         for g = 1:numg
%         	S(g,i) = in.det(g,i)*in.xsec((in.mt(i)-1)*numg+g, 1);
%         end
%     end
end

% -------------------------------------------------------------------------
% -------------------------------------- Discretizations
[mu,w] = S_1D(in.ord);
j = 0;
for i = 1:length(in.xfm)
    dx( (j+1):(j+in.xfm(i))   )  = (in.xcm(i+1)-in.xcm(i))/in.xfm(i);
    if nargin == 2
        for g=gbound(1):git:gbound(2)
            Q( (j+1):(j+in.xfm(i)), :, g)  = S(g,i);
        end
    else
        for g=gbound(1):git:gbound(2)
            for k = 1:ord
                Q( (j+1):(j+in.xfm(i)),k,g)  = qad((j+1):(j+in.xfm(i)),g);
            end
        end
    end
    mtt( (j+1):(j+in.xfm(i))   )  = in.mt(i);  % assign mat to each f mesh
    j = sum(in.xfm(1:i));
end
n = sum(in.xfm);

% -------------------------------------------------------------------------
% -------------------------------------- Matrix Pre-allocation
psi         =   zeros(n+1,ord,numg);        % Angular Flux
cur         =   zeros(n,numg);              % Current
phi_o       =   zeros(n,numg);              % Scalar Flux Old
phi         =   zeros(n,numg);              % Scalar Flux New
con1        =   zeros(n,ord,numg);          % Transport Coef.
con2        =   zeros(n,ord,numg);          % Transport Coef.
con3        =   zeros(n,ord,numg);          % Transport Coef.
q           =   zeros(n,ord,numg);          % Scattered Source

% -------------------------------------------------------------------------
% -------------------------------------- Transport Coefficients
for g = gbound(1):git:gbound(2)
    for k = 1:n
        m = mtt(k);
        con1(k,:,g) = in.xsec((m-1)*numg+g,1)*dx(k)./(2.0*mu(:));
        con2(k,:,g) = (1.0-sign(mu(:)').*con1(k,:,g))./...
            (1.0+sign(mu(:)').*con1(k,:,g));
        con3(k,:,g) = sign(mu(:)').*dx(k)./(mu(:)'.*...
            (1.0+sign(mu(:)').*con1(k,:,g)));
    end
end

q = Q;

% ----------------- Solution Algorithm ------------------------------------

for g = gbound(1):git:gbound(2)
    
    % ----------------- Convergence Parameters ----------------------------
    err  = 10;
    iter = 0;
    
    while err > in.maxerr && iter < in.maxit
        % -------------------------------------- Mu > 0
        for k = 1:1:n
            for j = bound(1,1):bound(1,2)
                psi(k+1,j,g) = psi(k,j,g)*con2(k,j,g)+q(k,j,g)*con3(k,j,g);
            end
        end
        % -------------------------------------- Mu < 0
        for k = n:-1:1
            for j = bound(2,1):bound(2,2)
                psi(k,j,g) = psi(k+1,j,g)*con2(k,j,g)+q(k,j,g)*con3(k,j,g);
            end
        end
        
        % ---- Scalar Flux
        for k = 1:n
            phi(k,g) = sum(w'.*0.25.*(psi(k,:,g)+psi(k+1,:,g)));
        end
        % ---- Updated Source Term
        for z   =   g:git:gbound(2) % only down scattering
            for kk = 1:n
                q(kk,:,z) = Q(kk,:,z); % reset
            end
            if adj == 0
                for k = 1:n
                    m = mtt(k);
                    for gg = gbound(1):git:gbound(2) % group gg to group z
                        q(k,:,z) = q(k,:,z) + ...
                            in.xsec((m-1)*numg+gg,2+z)*phi(k,gg);
                    end
                end
            else
                for k = 1:n
                    m = mtt(k);
                    for gg = gbound(1):git:gbound(2) % group gg to group z
                        q(k,:,z) = q(k,:,z) + ...
                            in.xsec((m-1)*numg+z,2+gg)*phi(k,gg);
                    end
                end
            end
        end
        % ---- Scalar Flux Error Between Iterations
        err = max(max(abs(phi-phi_o)./phi));
        % ---- Reset Scalar Flux
        phi_o = phi;
        % ---- Iteration Counter
        iter = iter + 1;
    end
    
end

% stuff for future plotting: phiPLOT is a step-wise phi, and phiAVG is the
% average value for a coarse mesh region (for comparison to MC values)


xa(1)=dx(1)/2+in.xcm(1);
for i = 2:length(dx)
    xa(i) = xa(i-1) + 0.5*(dx(i-1)+dx(i));
end
xe(1)=in.xcm(1);
for i = 2:length(dx)+1
    xe(i) = xe(i-1) + dx(i-1);
end

    
for g = 1:in.numg
    j=0;
    for i = 1:length(in.xfm)
        phiPLOT( (j+1):(j+in.xfm(i)), g)  = ...
            mean( phi( (j+1):(j+in.xfm(i)) , g) ) ;
        phiAVG(i,g)   = mean( phi( (j+1):(j+in.xfm(i)) , g) );
        j = sum(in.xfm(1:i));
    end
end

if sum(sum(phi<0)) > 0
    if (adj==1)
        disp('*** warning: negative adjoint')
    else
        disp('*** warning: negative forward')
    end
    disp('     make smaller mesh size or turn on')
    disp('     negative flux fixup')
end
if sum(sum(sum((psi<0)))) > 0
    if (adj==1)
        disp('*** warning: negative angular adjoint')
    else
        disp('*** warning: negative angular forward')
    end
    disp('     make smaller mesh size or turn on')
    disp('     negative flux fixup')
end

end