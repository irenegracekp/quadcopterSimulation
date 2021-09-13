start_time = 0;
end_time = 10;
dt = 0.05;
times = start_time:dt:end_time;

% Number of points in the simulation.
N = numel(times);

% Initial simulation state.
x = [0; 0; 10];
xdot = zeros(3, 1);


t = start_time;
i = 1;
m = 1; %kg
l = 0.1; %meter
J = [0.01 0 0; 
     0 0.05 0; 
     0 0 0.001];
position_initial = [0; 0; 10];
theta = zeros(3, 1);
Positiondot_initial = zeros(3, 1);

% Simulate some disturbance in the angular velocity.
% The magnitude of the deviation is in radians / second.
deviation = 100;
thetadot = zeros(3, 1);
theta = deg2rad(2 * deviation * rand(3, 1) - deviation);

U = zeros(4, 1);
ThetaRef = zeros(3, 1);
realTheta = 0;
while( i < 55)
    
    [thetadot,thetadotdot] = controller(theta, thetadot, ThetaRef) ;
    Moments = toMoments(thetadotdot,J);
    
    theta = theta + (thetadot*dt);
    thetadot = thetadot + (thetadotdot*dt);
    
    realTheta_1(i) = thetadot(1);
    realTheta_2(i) = thetadot(2);
    realTheta_3(i) = thetadot(2);
    %omega = rotation_itobody * thetadot ; 
    

    
    i = i+1;
end
t = tiledlayout(1,3); % Requires R2019b or later
nexttile
plot(realTheta_1)
nexttile
plot(realTheta_2)
nexttile
plot(realTheta_3)

t.Padding = 'compact';
t.TileSpacing = 'compact';


function [thetadot,thetadotdot] = controller(curTheta, curThetaDot, ThetaRef)
    
    kp = 30;
    kd = 0.1;
    
    e = ThetaRef - curTheta;
    ThetadotRef = kp*e ;
    edot = ThetadotRef - curThetaDot;
    
    thetadot = kp*e ;
    thetadotdot = kd*edot;
      
end

function Moments = toMoments(thetadotdot,J)
    
    Moments = thetadotdot'/J;
    Moments =  Moments'
        
end



