% Hoogh_2018 Initially used for Groundwater Flow Model
% Based on Hooghoudt-situation
% Modelling and Simulating 2018
% Author: Willem Bouten (2015)

% Changed to represent an Erosion Model (after Minasny etal (2001))
% Applied to an area in Luxemburg
% By: Marrit Leenstra (2018)

%%%%%%%%%%%%% INITIALISATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

load DEM_Lux.txt

% Control constants
PixX = 25; PixY = 25;           % Pixel length in X and Y-direction 	[m]
[ny,nx]=size(DEM_Lux); 			% Number of Pixels in X and Y-direction [-]
StartTime = 0;					% Start Time for simulation 			[year]
EndTime = 2000;			        % Time at which simulation Ends 		[year]
dt = 1;    					% calculation time step 				[year] 

% System constants after Minasny and McBratney (2001)
D = 0.08;                      % Erosive diffusitivity of material     [L2/T]
P0 = 1E-3;                      % Potential physical weathering rate    [m/year]
b = 1.5;                        % Empirical constant weathering         [/L]
W0 = 1.6E-3;                    % Potential chemical weathering rate    [m/year]
k1 = 0.8;                       % Rate constant for soil thickness      [/L]
Pr = 2.6;                       % Rockdensity                           [Mg/m3]
Ps = 1.6;                       % Soildensity                           [Mg/m3]

% Initialisation and boundary conditions
Time = StartTime;				% Initialisation of Time				[day]
h(1:ny,1:nx) = 2;               % Initial soil thickness                [m]
DEM_soil = DEM_Lux + h;         % Digital Elevation Model surfacelayer  [m]
DEM_rock = DEM_Lux;             % Digital Elevation Model bedrock	    [m]
FlowX(:,1:nx+1) = 0; 
FlowY(1:ny+1,:) = 0;            % no flow (= erosion)                   [L3 L-1 t-1]

% Initialisation soil balance variables
TotSoil = sum(h(:));
TotErosion = 0;

% Auxillary constants
x = [1:nx]*PixX-PixX/2;         % mid of x-grid                         [m]
y = [1:ny]*PixY-PixY/2;         % mid of y-grid                         [m]

% %%%%%%%%%%%% DYNAMIC LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%
while Time <= EndTime
 
% Calc. Flow (erosion) in x-direction
for j = 2:nx-1
    for i = 1:ny
        %DAX(i,j) = D * PixY * (h(i,j)+h(i,j-1))/2;
        FlowX(i,j) = D*((DEM_soil(i,j-1)-(2*DEM_soil(i,j))+DEM_soil(i,j+1))/ ...
            (PixX^2))*dt;
    end
end

% Calc. Flow (erosion) in y-direction
for j = 1:nx
    for i = 2:ny-1
        %DAY(i,j) = D * PixY * (h(i,j)+h(i-1,j))/2;
        FlowY(i,j) = D*((DEM_soil(i-1,j)-(2*DEM_soil(i,j))+DEM_soil(i+1,j))/ ...
            (PixY^2))*dt;
    end
end

% Avoid negative soiltickness
[FlowX, FlowY] = zerocheck((DEM_soil - DEM_rock), FlowX, FlowY);

% Calculate chemical weathering
for j = 1:nx
    for i = 1:ny
        if (h(i,j) > 0)
            W(i,j) = W0*(1-exp((-k1*h(i,j))))*dt;
        end
    end
end

% Calculate physical weathering
for j = 1:nx
    for i = 1:ny
        if h(i,j) > 0
            P(i,j) = -P0*exp(-b*h(i,j))*dt;
        end
    end
end

% Add rockdensity and soildensity to physical weathering rate
E(i,j) = -(Pr/Ps)*P(i,j);

% Calc h with forward integration
for j = 2:nx
    for i = 2:ny
        netFlow(i,j)= FlowX(i,j-1)-FlowX(i,j+1)+FlowY(i-1,j)-FlowY(i+1,j);
        h(i,j) = h(i,j) + netFlow(i,j) ;
    end
end

% Update DEM_soil (surface)
DEM_rock = DEM_rock - E(i,j);
DEM_soil = DEM_rock + E(i,j) + h - W(i,j);
h = DEM_soil-DEM_rock;
% DEM_soil = DEM_rock + h;





% % Calc. Flow in x-dir: Flow = -KA * dH/dx; (KA = K * average height of watercolumn * PixY)
% KAX = K * PixY * (H(:,1:nx-1)-DEM_Clay(:,1:nx-1) + H(:,2:nx)-DEM_Clay(:,2:nx))/2;     %	[m3/dag]
% FlowX(:,2:nx) = -1 * KAX .* (H(:,2:nx)-H(:,1:nx-1)) / PixX;                           %	[m3/day]
%  
% % Calc. Flow in y-dir: Flow = -KA * dH/dy; (KA = K * average height of watercolumn * PixX)
% KAY = K * (0.5*(H(1:ny-1,:)-DEM_Clay(1:ny-1,:) + H(2:ny,:)-DEM_Clay(2:ny,:))) * PixX;   %	[m3/dag]
% FlowY(2:ny,1:nx) = -1 * KAY .* (H(2:ny,:)-H(1:(ny-1),:)) / PixY;                        %	[m3/day]
% 
% % Calculate precipitation rate in Volume per gridcell;
% PrecRVol = PrecR * PixX * PixY;                                                         % 	[m3/day]
% 
% % Calculate new H values by forward integration
% NetFlow(1:ny,2:nx-1) = FlowX(1:ny,2:nx-1) - FlowX(1:ny,3:nx) + ...
%                              FlowY(1:ny,2:nx-1) - FlowY(2:ny+1,2:nx-1);                 %	[m3/day]
% H(1:ny,2:nx-1) = H(1:ny,2:nx-1) + ...
%                        (PrecRVol(1:ny,2:nx-1)+ NetFlow(1:ny,2:nx-1))* dt/ ...
% 					   (PorVol*PixX*PixY);                                              %	[m]
 Time = Time + dt;
 
 if mod(Time,100)==0
        figure(10); 
        set(gcf,'Position', [0 0 1600 1600]);
        title ([num2str(Time)])
        %subplot(1,2,1),
        imagesc(h, [0 10]); colorbar; title(['Soil Thickness,  time = ' ,num2str(Time)]);

        %subplot(1,2,2),
        %imagesc(e,[150 450]); colorbar; title('Dynamic DEM');

        drawnow;
 
%   if mod(Time,1000)< dt,				% plot every 100 years
     
%      figure(5);
%      plot(h(:,100));
%      hold on;
%    figure(5); contourf(x,y,h,20); colorbar;
%    title('From above'); xlabel('distance [m]'); ylabel('distance [m]');
%    figure(4); 
%    plot(DEM_soil(:,100), 'r'); 
%    hold on
%    plot(DEM_rock(:,100), 'k');
%    plot(h(:,100), 'g')
%    legend('surface', 'bedrock', 'soil-rock')
%    title('Cross section'); 
%    xlabel(' distance [m]'); ylabel('hight [m]');
    drawnow
  end;
% 
 end %while Time <= EndTime
% %%%%%%%%%%%% END DYNAMIC LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%

% calc. soil balance
TotErosion = TotErosion + sum(FlowX(:)) + sum(FlowY(:));
TotSoil = TotSoil - TotErosion;

%figure(1)
%contourf(DEM_soil); colorbar; xlabel('Nr gridcell x-direction'); ylabel('Nr gridcell y-direction')

figure(1)
contourf(h);
colorbar;
title('From above'); 
xlabel('distance [m]'); ylabel('distance [m]');

figure(2)
%mesh(rot90(DEM_soil,3)); 
%hold on
mesh(h); colorbar;

figure(3)
plot(DEM_soil(:,100), 'r'); 
hold on
plot(DEM_rock(:,100), 'k');
plot(h(:,100)*10, 'g');
legend('surface', 'bedrock');
hold off