function output = mcslab1d(in,srcout,vrout)
% function output = mcslab1d(in,vrout)
%   This function performs multi-group Monte Carlo particle transport in a
%   1-D Cartesian slab.
%
% J. Roberts 4/13/2010

	tic
	
	% initialize flux estimators
	trakest = zeros(2, in.numslabs, in.numg);
	collest = zeros(2, in.numslabs, in.numg);
    % monte carlo particle population estimator
    if in.mcparts == 1
        mcpest = zeros(2, in.numslabs, in.numg);
    end
    
    ww = in.stcadis+in.fwcadis+in.coopers; % >0 to check ww's

	minw = 1;
    
    %%%%%%%%%%%%%%%%% BEGIN HISTORIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for n = 1:in.N % histories

        % Call source() to get my starting location and direction.  Because
        % source biasing may yield several particles, I've got to have a
        % bank of neutrons associated with history "n".  This bank will be
        % useful for geometry splitting and weight windows, too.
        [x,mu,g,ps,wt] = source(in,srcout,vrout);
        bank(1:ps,:) = [ x mu g wt ];
        
        % reset (or initialize) temporary tallies
        tmptrak = zeros(1,in.numslabs,in.numg);
        tmpcoll = zeros(1,in.numslabs,in.numg);
        tmpmcp  = zeros(1,in.numslabs,in.numg);
        
        while length(bank(:,1)) >= 1 % go through entire bank
           % if length(bank(:,1))>10000
           %     disp('whoa!')
           %     return
           % end
            x  = bank(1,1);  % my start location
            mu = bank(1,2);  % my start direction
            g  = bank(1,3);  % my energy group
            w  = bank(1,4);  % my weight
            alive = 1;       % i am alive
            coll = 0;        % i am not exiting a collision
            mycell = getmycell(x,in);  % get my cell
            
            while alive == 1;
                
                % post collision stuff
                if coll == 1
                    % i need a new direction in life
                    mu = 2*rand-1; % we're isotropic
                    % i also need a new energy 
                    % sct gives me the vector of g->g' sigs where I am g
                    sct = in.xsec( in.numg*(in.mt(mycell)-1)+g, 3:end);
                    sct = sct/sum(sct); % normalize for prob.
                    ksi = rand;
                    for i = 1:length(sct)
                        if ksi < sum(sct(1:i))
                            g = i;
                            break
                        end
                   end           
                end
                
                % check weight window if needed
                if ww>0
                   [w,bank,alive] = wwcheck(x,mu,g,w,bank,vrout);
                   if alive == 0
                       break
                   end
                end   
                
                % determine my neighbor and how far to her
                [neighbor,d2neighbor] = getmyneighbor(mycell,x,mu,in);
                
                % sample how many mfp's I will go
                ksi = rand;
                mfps = -log(ksi);
                
                % determine my total cross-section
                SigT = in.xsec( in.numg*(in.mt(mycell)-1)+g, 1);
                
                % determine how far along x-axis I go
                dist = mfps*mu/SigT;
                
                % determine whether I reach surface
                if ( abs(dist) > d2neighbor) % then I pass to the next cell
                    tmptrak(1,mycell,g)=tmptrak(1,mycell,g)+...
                        w*abs(d2neighbor/mu);
                    if  (neighbor==0)
                        % I leaked out of left
                        alive=0;
                    elseif (neighbor==in.numslabs+1)
                        % I leaked out of right
                        alive=0;
                    else % i simply get to the boundary
                        % and i may or may not be split or rouletted
                        x = x+sign(mu)*d2neighbor;
                        if in.geosplt == 1 && neighbor > 0 ...
                                           && neighbor < in.numslabs+1     
                          [w,bank,alive] = ...
                              splitme(neighbor,mycell,x,mu,g,w,bank,vrout);
                        end
                        mycell = neighbor;
                        coll = 0;
                    end
                else
                    % otherwise, I collide; no matter what collision we
                    % have, tally my incoming weight
                    tmpcoll(1,mycell,g) = tmpcoll(1,mycell,g)+w;
                    tmptrak(1,mycell,g) = tmptrak(1,mycell,g)+w*mfps/SigT;
                    if in.mcparts == 1
                    	tmpmcp(1,mycell,g)  = tmpmcp(1,mycell,g)+1; 
                    end
                    
                    % update my location
                    x = x+dist;
                    if in.impcapt == 1
                        % reduce my weight
                        w = w*(1-in.xsec(in.numg*(in.mt(mycell)-1)+g,2)...
                            /SigT);
                        coll = 1;
                        % check to see if weight is below cutoff
                        if w < in.wcut
                            if rand < w/in.wavg
                                w = in.wavg;
                            else
                                alive = 0;
                                break;
                            end
                        end
                    else
                        if rand < in.xsec(in.numg*(in.mt(mycell)-1)+g,2)...
                                /SigT;
                            % i am absorbed
                            alive = 0;
                        else
                            % i am scattered;
                            coll = 1;
                        end
                    end
                end

            end
            % end while alive
            
            if length(bank(:,1)) > 1 % reduce the bank
                bank = bank(2:end,:);
            else % else, we'll go to the next history
                break
            end
            
        end
        % end while bank
        
        % update tallies
        for i = 1:in.numg
            collest(1,:,i) = collest(1,:,i) + tmpcoll(1,:,i);
            collest(2,:,i) = collest(2,:,i) + tmpcoll(1,:,i).^2;
            trakest(1,:,i) = trakest(1,:,i) + tmptrak(1,:,i);
            trakest(2,:,i) = trakest(2,:,i) + tmptrak(1,:,i).^2;
            if in.mcparts==1
                mcpest(1,:,i) = mcpest(1,:,i) + tmpmcp(1,:,i);
                mcpest(2,:,i) = mcpest(2,:,i) + tmpmcp(1,:,i).^2;
            end
        end

    end
    %%%%%%%%%%%%%%%%% END OF HISTORIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	%disp(['minw=',num2str(minw)])  
	t2 = toc;
	output = struct(  'trakest',    trakest, ...     % track length est
                      'collest',    collest, ...     % collision est
                            't',          t2 ...     % finish time
	);    
    if in.mcparts==1
       output.mcpest = mcpest; 
    end
end
% end function playgame

function [x,mu,g,ps,wt] = source(in,srcout,vrout)
	% In this function, the source distribution is assigned.
 
    % Analog Source
    %   Here, we use the pre-generated source distribution vectors
    %   from srcdriver.  Px is the discrete pdf if a particle is born
    %   (uniformly) within a given slab cell.  Once that cell is
    %   selected, x is sampled uniformly within the cell.  Then
    %   the group is sampled within that cell using 
    %      p(group|cell) = p(group,cell)/p(cell) = p(:,i)/Px(i)
    if (in.srcbias + in.stcadis + in.fwcadis + in.coopers) == 0
        
        g=-1; x=-1;
        ksi1 = rand;
        for i = 1:srcout.numsrc
            if ksi1 < sum(srcout.Px(1:i))
                x = rand*(srcout.sx(i,2)-srcout.sx(i,1))+srcout.sx(i,1);
                break
            end
        end
        ksi2 = rand;
        for j = 1:in.numg
            if ksi2 < sum(srcout.P(1:j,i))/srcout.Px(i)
                g = j;
                break
            end
        end
        if g==-1 || x == -1
            disp('kaka')
        end
        mu = 2*rand-1;
        wt = 1;
        ps = 1;

    elseif ( in.stcadis + in.fwcadis + in.coopers ) > 0
        
        g=-1; x=-1;
        ksi1 = rand;
        for i = 1:vrout.numsrc
            if ksi1 < sum(vrout.Pxf(1:i))
                x = rand*(vrout.sxf(i,2)-vrout.sxf(i,1))+vrout.sxf(i,1);
                break
            end
        end
        if in.numg > 1
            ksi2 = rand;
            for j = 1:in.numg
                if ksi2 < sum(vrout.Pf(1:j,i))/vrout.Pxf(i)
                    g = j;
                    break
                end
            end
        else
            g = 1;
        end
        if g == -1 || x == -1
            disp('source issue!')
        end
        mu = 2*rand-1;
        wt = vrout.weight(g,i);
        ps = 1;
    end
        
    % MANUAL SOURCE BIASING 
    %   feel free to add me
end
% end function source

function [neighbor, d2neighbor] = getmyneighbor(mycell,x,mu,in)
    % mycell is my current cell, and neighbor is the one i'd enter
    if mu > 0
        d2neighbor = (in.xcm(mycell+1)-x);
        neighbor=mycell+1;
    else
        d2neighbor = (x-in.xcm(mycell));
        neighbor=mycell-1;
    end
end
% end function getmyneighbor

function [mycell] = getmycell(x,in)
    for mycell = 1:in.numslabs
        if x < in.xcm(mycell+1)
            break
        end
    end
end
% end function mycell

function [w,bank,alive] = wwcheck(x,mu,g,w,bank,vrout)
    alive = 1;
    % first find my wL
% linear search:  (r
%     for i = 1:length(vrout.xww(:,1))
%         if x < vrout.xww(i,2)
%             wL = vrout.wL(i,g);
%             wU = vrout.cU*wL;
%             break
%         end
%     end
    i = binsearch(x,vrout.xww(:,2),mu); %mu used to correct index at bound
    wL = vrout.wL(i,g); wU = vrout.cU*wL;
    if w > wU
%         r = w/wU; % NEW
%         if rand < r-floor(r)
%             n = ceil(r);
%         else
%             n = floor(r);
%         end
%         w = w/n;
%         if n > 1
%           bank(end+1:end+n-1,:) = ...
%             [ ones(n-1,1)*x ones(n-1,1)*mu ones(n-1,1)*g ones(n-1,1)*w ];
%         end
        n = floor(1+w/wU); % OLD
        if n > 1
            w = w/n;
            % bank n-1 particles
            bank(end+1:end+n-1,:) = ...
             [ ones(n-1,1)*x ones(n-1,1)*mu ones(n-1,1)*g ones(n-1,1)*w ];
        end
    elseif w < wL
        P = 2*w/(wU+wL); % i.e. the average w is in the center
        if rand < P
            w = w/P;
        else
            alive = 0;
        end
    end
end
% end function wwcheck

function index = binsearch(x,F,mu)
    % simple binary search for use in source sampling and wwcheck
    L = 1;  R = length(F);
    while(1)
        M = round( (L+R)/2 );
        if x > F(M)
            L = M+1;
        elseif x < F(M)
            R = M-1;
        else
            index=M+(mu>0);
            break
        end
        if ( L>R)  % end of search; determine which side we end on
            if F(M) > x
                index=M;
            else
                index=M+1;
            end
            break
        end
    end
end
% end function binsearch

function [w,bank,alive] = splitme(neighbor,mycell,x,mu,g,w,bank,vrout)
    alive = 1;
    r = vrout.imp(neighbor,g)/vrout.imp(mycell,g);
    if (r > 1.0) % split
        if rand < r-floor(r)
            n = ceil(r);
        else
            n = floor(r);
        end
        w = w/r;
        if n > 1
          bank(end+1:end+n-1,:) = ...
            [ ones(n-1,1)*x ones(n-1,1)*mu ones(n-1,1)*g ones(n-1,1)*w ];
        end
    else % russian roulette
        if rand < r
            w = w/r;
        else
            alive = 0;
        end
    end
end
% end function splitme

