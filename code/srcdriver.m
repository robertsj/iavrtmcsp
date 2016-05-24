function srcout = srcdriver(in)
% function srcout = srcdriver(in)
%   This function produces the unbiased source density
%   for use in the simulation.  Its basic output is the
%   cell/energy discrete pdf.

% source biasing
    %  Here, we need to look at the volumetric source distribution
    %  as well as the fine mesh adjoint
    if in.numg>1
        srcloc = find(sum(in.src)>0); % picks out cell numbers w/ src
    else
        srcloc = find(in.src>0);
    end
    denom = 0;
    
    for i = 1:length(srcloc)
        V(i) = in.xcm( srcloc(i)+1)-in.xcm( srcloc(i) );
        sx(i,1)=in.xcm( srcloc(i) ); % lower bound
        sx(i,2)=in.xcm( srcloc(i)+1); % upper bound
        denom = denom + V(i)*sum( in.src(:,srcloc(i)));  
    end
    for i = 1:length(srcloc)
        P(1:in.numg,i)=V(i)*in.src(1:in.numg, srcloc(i))/denom;    
    end
    % P gives a cell/energy-dependent source density
    %   To use this in an unbiased source, so the following:
    %       Px = sum(P) --> cell dependent probability
    %       for i=1,numsrccell, if rand < sum(Px(1:i)), found cell
    %        x = rand*(dx(i,2)-dx(i,1))+dx(i,1)
    %       for j=1,numg, if rand < sum(Px(i,j)), found group
    %        g = j;
    if in.numg > 1
        Px=sum(P);
    else
        Px = P;
    end
    numsrc=length(srcloc);
    
    srcout   =   struct( 'srcloc', srcloc, ...
                         'numsrc', numsrc, ...
                         'sx',     sx, ...
                         'P',      P, ...
                         'Px',     Px, ...
                         'nrmfct', denom);
end