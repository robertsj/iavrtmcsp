function finout = out1d(in,out,snout,srcout,vrout)
% function out1d(in,out,snout)
%   This function manipulates the tallies from the simulation into the
%   fluxes or reaction rates required.  It computes figure-of-merits.  It
%   plots the Monte Carlo results and the (forward) discrete ordinates
%   results.

	disp(['  Monte Carlo elapsed time: ', num2str(out.t)])
    disp(['   Forward SN elapsed time: ', num2str(snout.snt)])
    if ( in.geosplt+in.stcadis+in.fwcadis+in.coopers > 0 )
        if in.fwcadis == 1
            vrout.t=vrout.t+snout.snt;
        end
        if in.coopers == 1
            vrout.t=vrout.t+snout.snt;
        end
        disp(['    VR Param. elapsed time: ', num2str(vrout.t)])
    end
    disp(' Variance Reduction: ')
    if in.impcapt == 1, disp(' * implicit capture'), end
    if in.srcbias == 1, disp(' * manual source biasing'), end
    if in.geosplt == 1, disp(' * adjoint-based geometry splitting'), end
    if in.stcadis == 1, disp(' * standard cadis'), end
    if in.fwcadis == 1, disp(' * forward-weighted cadis'), end
    if in.coopers == 1, disp(' * pseudo coopers'), end
    
    [phi1,re1]   = phi_tl(in,out);
    [phi2,re2]   = phi_cd(in,out);

%     if in.mcparts == 1
%         out.collest = out.mcpest;
%         [phiMC,reMC] = phi_cd(in,out);
%     end
    
    if in.mcparts == 1
    	[phiMC,reMC] = phi_mc(in,out);
    end
        
    finout.t=out.t+vrout.t;
    % now print the fluxes all nice
    for g = 1:in.numg
    disp([' group ',num2str(g)])
    disp(' reg   |     phi1       RE         FOM    |     phi1       RE         FOM   | ')
    disp('-----------------------------------------------------------------------------')
    for i = 1:in.numslabs
        if in.det(g,i) > 0 % point an arrow at a detector tally
        fprintf(1,'%4i   |  %10.4e %10.4e %5.3e | %10.4e %10.4e %5.3e |<-- \n', ...
               i, phi1(i,g), re1(i,g), 1/(re1(i,g)^2*(out.t+vrout.t)), ...
                  phi2(i,g), re2(i,g), 1/(re2(i,g)^2*(out.t+vrout.t))  ); 
        else
        fprintf(1,'%4i   |  %10.4e %10.4e %5.3e | %10.4e %10.4e %5.3e | \n', ...
               i, phi1(i,g), re1(i,g), 1/(re1(i,g)^2*(out.t+vrout.t)), ...
                  phi2(i,g), re2(i,g), 1/(re2(i,g)^2*(out.t+vrout.t))  );
        end
    end
    disp('-----------------------------------------------------------------------------')

    end
    
    xx = 0.5*(in.xcm(1:end-1)+in.xcm(2:end));
    figure(1)
    hold on
    for i = 1:in.numg
        fig = errorbar(xx,phi2(:,i),re2(:,i).*phi2(:,i),'o');
        set(fig,'Color',plotcolor(i),'MarkerEdgeColor','k',...
                            'MarkerFaceColor',plotcolor(i),...
                            'MarkerSize',6);
        lab(i,:)=(['\phi_',num2str(i)]);    
    end
    for i = 1:in.numg 
        % Sn plot (normalized)
        plot(snout.x, snout.phiFplot(:,i)/srcout.nrmfct, ...
            'Color',plotcolor(i),'LineStyle','--','LineWidth',2) 
    end
    %axis([in.xcm(1) in.xcm(end) 0 1.25*max(max(phi2))])
    title('Normalized Flux (S_N approx. dashed)')
    xlabel('x [cm]'), ylabel('\phi(x,E_g) [n/cm^2-sp]')
    grid on
    legend(lab,0)
    
    % final output
    finout.phi1 = phi1; finout.re1 = re1;
    finout.phi2 = phi2; finout.re2 = re2;
    finout.fom1 = 1./(re1.^2*finout.t);
    finout.fom2 = 1./(re2.^2*finout.t);
    if in.mcparts == 1
        finout.phiMC = phiMC; finout.reMC = reMC;
    end
   
    figure(2)
    title('Track Length Estimated Flux -- Relative Error')
    xlabel('x [cm]'), ylabel('relative error')
    hold on
    if in.impcapt==1
        plot(xx,re2,'Color',plotcolor(1),'LineWidth',2)
    elseif in.geosplt == 1
        plot(xx,re2,'Color',plotcolor(2),'LineWidth',2)
    elseif in.stcadis == 1
        plot(xx,re2,'Color',plotcolor(3),'LineWidth',2)       
    elseif in.fwcadis == 1
        plot(xx,re2,'Color',plotcolor(5),'LineWidth',2)
    elseif in.coopers == 1
        plot(xx,re2,'Color',plotcolor(6),'LineWidth',2)
    end
    grid on
    
    if in.mcparts==1
        figure(3)
        title('Monte Carlo Particle Density')
        xlabel('x [cm]'), ylabel('m(x,E_g) [mcp/cm^3-sp]')
        hold on
        for i = 1:in.numg
            fig = errorbar(xx,phiMC(:,i),reMC(:,i).*phiMC(:,i),'o-');
            set(fig,'Color',plotcolor(i),'MarkerEdgeColor','k',...
                'MarkerFaceColor',plotcolor(i),...
                'MarkerSize',6);
            lab(i,:)=(['\phi_',num2str(i)]);
        end
    end
    grid on
    
end

function [phi,re] = phi_tl(in,out)
    % compute track length estimate of flux
    S1 = zeros(in.numslabs,in.numg);
    S2 = S1; phi = S1; re = S1;
    V  = in.xcm(2:end)-in.xcm(1:end-1);
    for g = 1:in.numg
        S1(:,g)   = out.trakest(1,:,g)';
        phi(:,g) = S1(:,g)./(in.N*V)';
        S2(:,g)   = out.trakest(2,:,g)';
        re(:,g)   = sqrt(  1/(in.N-1)' * (S2(:,g)./(in.N*V.^2)' - ...
                     (phi(:,g)).^2) )./phi(:,g);
    end 
end

function [phi,re] = phi_cd(in,out)
    % compute collision density estimate of flux
    % compute track length estimate of flux
    S1 = zeros(in.numslabs,in.numg);
    S2 = S1; phi = S1; re = S1;
    V  = in.xcm(2:end)-in.xcm(1:end-1);    
    for g = 1:in.numg
        S1(:,g) = out.collest(1,:,g)';
        sigT    = zeros(1,in.numslabs);
        for i = 1:in.numslabs
            sigT(i) = in.xsec( in.numg*(in.mt(i)-1)+g, 1);
        end
        sig_v    = sigT.*V;
        phi(:,g) = S1(:,g)./(in.N*sig_v)';
        S2(:,g)  = out.collest(2,:,g)';
        re(:,g)  = sqrt(  1/(in.N-1)' * ( S2(:,g)./(in.N*(sig_v).^2)'- ...
                    phi(:,g).^2 ) )./phi(:,g);
    end
end

function [phi,re] = phi_mc(in,out)
    % compute collision density estimate of flux
    % compute track length estimate of flux
    S1 = zeros(in.numslabs,in.numg);
    S2 = S1; phi = S1; re = S1;
    V  = in.xcm(2:end)-in.xcm(1:end-1);    
    for g = 1:in.numg
        S1(:,g) = out.mcpest(1,:,g)';
        sig_v    = V;
        phi(:,g) = S1(:,g)./(in.N*sig_v)';
        S2(:,g)  = out.mcpest(2,:,g)';
        re(:,g)  = sqrt(  1/(in.N-1)' * ( S2(:,g)./(in.N*(sig_v).^2)'- ...
                    phi(:,g).^2 ) )./phi(:,g);
    end
end


function color = plotcolor(g)
    % this function is a hard-coded color map for the different
    % flux groups.  I've accounted for up to 8 groups.
    % set(ur,'Color',[1 0.7 0.2],'LineWidth',2);
    switch g
        case 1
            color = [0.0 0.0 1.0]; % blue
        case 2
            color = [0.0 0.8 0.2]; % nice green
        case 3
            color = [1.0 0.0 0.0]; % red
        case 4
            color = [0.4 0.0 0.6]; % purple
        case 5
            color = [0.9 0.4 0.0]; % orange
        case 6
            color = [0.5 0.2 0.0]; % brown
        case 7
            color = [0.0 0.8 0.6]; % turquoise
        case 8
            color = [0.7 0.6 0.0]; % gold
        otherwise
            
    end
end