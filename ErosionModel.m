% Hoogh_2018 Initially used for Groundwater Flow Modele
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
EndTime = 400;			        % Time at which simulation Ends 		[year]
dt = 10;    					% calculation time step 				[year] 

% System constrants after Minasny and McBratney (2001)
D = 0.008;                      % Erosive diffusitivity of material     [L2/T]
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


% Auxillary constants
x = [1:nx]*PixX-PixX/2;         % mid of x-grid                         [m]
y = [1:ny]*PixY-PixY/2;         % mid of y-grid                         [m]

% %%%%%%%%%%%% DYNAMIC LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% while Time <= EndTime
% 
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
% Time = Time + dt;
% if mod(Time,10)< dt,				% plot every 10 days
%    subplot(2,1,1); contourf(x,y,H,20); 
%    title('From above'); xlabel('distance [m]'); ylabel('distance [m]');
%    subplot(2,1,2); plot(x,H(2,:)); ylim([0 1]); hold on 
%    title('Cross section'); 
%    xlabel(' distance [m]'); ylabel('Groundwater level [m]');
%    drawnow
% end;
% 
% end %while Time <= EndTime
% %%%%%%%%%%%% END DYNAMIC LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure(1)
%contourf(DEM_soil); colorbar; xlabel('Nr gridcell x-direction'); ylabel('Nr gridcell y-direction')

figure(2)
mesh(rot90(DEM_soil,3)); colorbar;

figure(3)
plot(DEM_soil(:,100), 'r'); 
hold on
plot(DEM_rock(:,100), 'k');
legend('surface', 'bedrock');
hold off