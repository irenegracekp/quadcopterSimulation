start_time = 0;
end_time = 10;
dt = 0.05;
times = start_time:dt:end_time;
vizRate = rateControl(2/dt);

t = start_time;
i = 1;
J = [0.01 0 0; 
     0 0.05 0; 
     0 0 0.001];
position = [0; 0; 10];
positiondot = zeros(3, 1);

% Simulate some disturbance in the angular velocity.
% The magnitude of the deviation is in radians / second.
deviation = 100;
thetadot = deg2rad(2 * deviation * rand(3, 1) - deviation) ;
theta = deg2rad(2 * deviation * rand(3, 1) - deviation);
omega = rotation_itobody(theta) * thetadot ; 

ThetaRef = zeros(3, 1);
realTheta = 0;

position_ref = [20; 1; 10];
positiondot_ref = zeros(3, 1);
j = 1;
k = 1;
while( j < 50)
    if (j<20)
        position_ref = [20; 1; 10];
        [thrust, ThetaRef, positiondot, positiondotdot] = position_controller(position, positiondot, position_ref, thetadot) ;
    elseif (j<40)
         position_ref = [10; 5; 20];
        [thrust, ThetaRef, positiondot, positiondotdot] = position_controller(position, positiondot, position_ref, thetadot) ;   
          
    else 
        position_ref = [30; 2; 10];
        [thrust, ThetaRef, positiondot, positiondotdot] = position_controller(position, positiondot, position_ref, thetadot) ;   
         
     end
        i = 1;
        ThetaRef(3) = 0;
        ThetaRef = ThetaRef';
    while(i < 10)
    
        [thetadot,thetadotdot] = controller(theta, thetadot, ThetaRef) ;
        Moments = toMoments(thetadotdot,J);

        omega = rotation_itobody(theta) * thetadot ; 
        
        thetadot = thetadot + (thetadotdot*dt);
        position = position + (positiondot*dt);
        positiondot = positiondot + (positiondotdot*dt);
        
        realTheta(k) = theta(1);
        realposition(k) = position(1) ; 

        i = i+1;
        theta = theta + (thetadot*dt);
        k = k + 1;
        
        plotTrVec = position;
        plot3(position(1),position(2),position(3));
        %plotRot = axang2quat([0 0 1 theta(3); 0 1 0 theta(2); 1 0 0 theta(1)]);
        plotRot = rotm2quat([1 sin(theta(1))*tan(theta(2)) cos(theta(1))*tan(theta(2)); 0 cos(theta(1)) -sin(theta(1)); 0 sin(theta(1))*sec(theta(2)) cos(theta(1))*sec(theta(2))]);
        plotTransforms(plotTrVec', plotRot, "MeshFilePath", "multirotor.stl", "meshcolor", [0 0 1],"View","3D", "FrameSize", 2);
        light;

        waitfor(vizRate);
        
    end
    j = j+1;
end


function [thetadot,thetadotdot] = controller(theta, thetadot, ThetaRef)
    
    kp = 30;
    kd = 0.1;
    
    e = ThetaRef - theta;
    ThetadotRef = kp*e ;
    edot = ThetadotRef - thetadot;
    
    thetadot = kp*e ;
    thetadotdot = kd*edot;
      
end

function Moments = toMoments(thetadotdot,J)
    
    Moments = thetadotdot'/J;
    Moments =  Moments' ;
        
end

function R = rotation_itobody(theta)
    
    R = [1 sin(theta(1))*tan(theta(2)) cos(theta(1))*tan(theta(2)); 0 cos(theta(1)) -sin(theta(1)); 0 sin(theta(1))*sec(theta(2)) cos(theta(1))*sec(theta(2))]; 

end

function [thrust, ThetaRef, positiondot, positiondotdot] = position_controller(position, positiondot, position_ref, omega) 

    kp = 3;
    kd = 0.1;
    
    e = position_ref - position;
    positiondotRef = kp*e ;
    edot = positiondotRef - positiondot;
    
    positiondot = kp*e ;
    positiondotdot = kd*edot;
    
    thrust = findThrust(positiondot, positiondotdot, omega);
    ThetaRef = findtheta(thrust) ; 
    
end

function thrust = findThrust(positiondot, positiondotdot, omega)
    
    m = 2;
    W = [0 -omega(3) omega(2);
        omega(3) 0 -omega(1);
        -omega(2) omega(1) 0];
    g = 9.81 ; 
    gravity = [0 0 m*g]';
    
    T = m*positiondotdot + m*W*positiondot - gravity ;
    thrust = T;
    
end

function ThetaRef = findtheta(T) 

    ThetaRef(1) = atan(T(2)/ T(3)) ;
    ThetaRef(2) = atan(-T(1)/(T(2)*sin(ThetaRef(1)) + T(3)*cos(ThetaRef(1))));
    
end