% Hoogh_2018
% Groundwater Flow Module
% Based on Hooghoudt-situation
% Modelling and Simulating 2018
% Author: Willem Bouten (2015)

%%%%%%%%%%%%% INITIALISATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

% Control constants
PixX = 0.25; PixY = 5.0; 		% Pixel length in X and Y-direction 	[m]
nx = 102; ny = 10;  			% Number of Pixels in X and Y-direction [-]
StartTime = 0;					% Start Time for simulation 			[day]
EndTime = 40;			        % Time at which simulation Ends 		[day]
dt = 0.005; 					% calculation time step 				[day] 

% System constants
DEM(1:ny,1:nx) = 0.75;          % Digital Elevation Model 			    [m]
TSand = 1.75;					% Thickness Sandlayer beneath DEM       [m]
K = 1; 							% hydraulic Conductivity sand			[m/day]
PorVol = 0.35; 				    % Pore Volume sand						[-]
DEM_Clay = DEM-TSand;           % Digital Elevation Model Clay          [m]

% Initialisation and boundary conditions
Time = StartTime;				% Initialisation of Time				[day]
H(1:ny,1:nx) = 0;               % Initial hydraulic Head (H) value 		[m]
% H(1:ny,1) = 0;                 % height in left ditch           	    [m]
% H(1:ny,nx) = 0;                % height in right ditch                 [m]
FlowY(1,:) = 0; FlowY(ny+1,:) = 0;% no flow in/out in Y-direction       [m3/day]
PrecR(1:ny,1:nx) = 0.01;        % field of Precipitation Rate 	        [m/day]

% Auxillary constants
x = [1:nx]*PixX-PixX/2;         % x-grid                                [m]
y = [1:ny]*PixY-PixY/2;         % y-grid                                [m]

%%%%%%%%%%%% DYNAMIC LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%

while Time <= EndTime

% Calc. Flow in x-dir: Flow = -KA * dH/dx; (KA = K * average height of watercolumn * PixY)
KAX = K * PixY * (H(:,1:nx-1)-DEM_Clay(:,1:nx-1) + H(:,2:nx)-DEM_Clay(:,2:nx))/2;     %	[m3/dag]
FlowX(:,2:nx) = -1 * KAX .* (H(:,2:nx)-H(:,1:nx-1)) / PixX;                           %	[m3/day]
 
% Calc. Flow in y-dir: Flow = -KA * dH/dy; (KA = K * average height of watercolumn * PixX)
KAY = K * (0.5*(H(1:ny-1,:)-DEM_Clay(1:ny-1,:) + H(2:ny,:)-DEM_Clay(2:ny,:))) * PixX;   %	[m3/dag]
FlowY(2:ny,1:nx) = -1 * KAY .* (H(2:ny,:)-H(1:(ny-1),:)) / PixY;                        %	[m3/day]

% Calculate precipitation rate in Volume per gridcell;
PrecRVol = PrecR * PixX * PixY;                                                         % 	[m3/day]

% Calculate new H values by forward integration
NetFlow(1:ny,2:nx-1) = FlowX(1:ny,2:nx-1) - FlowX(1:ny,3:nx) + ...
                             FlowY(1:ny,2:nx-1) - FlowY(2:ny+1,2:nx-1);                 %	[m3/day]
H(1:ny,2:nx-1) = H(1:ny,2:nx-1) + ...
                       (PrecRVol(1:ny,2:nx-1)+ NetFlow(1:ny,2:nx-1))* dt/ ...
					   (PorVol*PixX*PixY);                                              %	[m]
Time = Time + dt;
if mod(Time,10)< dt,				% plot every 10 days
   subplot(2,1,1); contourf(x,y,H,20); 
   title('From above'); xlabel('distance [m]'); ylabel('distance [m]');
   subplot(2,1,2); plot(x,H(2,:)); ylim([0 1]); hold on 
   title('Cross section'); 
   xlabel(' distance [m]'); ylabel('Groundwater level [m]');
   drawnow
end;

end %while Time <= EndTime
%%%%%%%%%%%% END DYNAMIC LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
contourf(H,[0:0.025:1]); colorbar; xlabel('Nr gridcell x-direction'); ylabel('Nr gridcell y-direction')