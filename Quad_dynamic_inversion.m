
pos = [0; 0; 10] ;
pos_dot = zeros(3, 1);

desired_pos = [20; 5; 10];
desired_pos_dot = zeros(3, 1);
time = 0;
time2 = 0;
vizRate = rateControl(2/0.05);
j = 1;
k = 1;
ttata = [0; 0; 0] ;
while (j < 50)
    
    if (j<20)
        desired_pos = [20; 1; 10];
    elseif (j<40)
        desired_pos = [10; 5; 10];
    else 
        desired_pos = [30; 2; 10];
    end

initial_e_pos = pos - desired_pos;
initial_e_dposdt = pos_dot - desired_pos_dot;

t_pos = 0:0.1:0.5;
[t_pos,x_error,y_error,z_error] = controller_pos(t_pos, initial_e_pos, initial_e_dposdt) ;

[Thrust, phi_desire,theta_desire ,x_dot,y_dot] = dyn_inv_pos( x_error,y_error,z_error, desired_pos, desired_pos_dot) ;
time = time + 1;
theta_desired = [phi_desire;theta_desire;0];
i = 1;
    while (i<10)
        
        e = ttata - theta_desired;

        theta_desired_dot = [0;0;0];
        theta_dot = [0;0;0];
        e_dot = theta_dot - theta_desired_dot;

        t = 0:1:1; 

        J = [0.01, 0, 0; 
             0, 0.05, 0; 
             0, 0, 0.001];

        initial_x = e;
        initial_dxdt = e_dot;

        [t,phi,theta,psi] = controller_att( t, initial_x, initial_dxdt,theta_desired,theta_desired_dot) ;

        moments = dyn_inv_att(t, phi, theta, psi, theta_desired, theta_desired_dot,J) ;
        time = time + 1;
        ttata = [phi(end,1);theta(end,1);psi(end,1)];
        pos(1) = pos(1) + x_dot;
        pos(2) = pos(2) + y_dot;

        realTheta(k) = theta(end,1);
        realposition(k) = pos(1) ;
        k = k + 1;
        time2 = time;
        i = i+1;
        
        
        plotTrVec = pos;
        %plotRot = axang2quat([0 0 1 theta(3); 0 1 0 theta(2); 1 0 0 theta(1)]);
        plotRot = rotm2quat([1 sin(phi(end))*tan(theta(end)) cos(phi(end))*tan(theta(end)); 0 cos(phi(end)) -sin(phi(end)); 0 sin(phi(end))*sec(theta(end)) cos(phi(end))*sec(theta(end))]);
        plotTransforms(plotTrVec', plotRot, "MeshFilePath", "multirotor.stl","View","3D", "FrameSize", 2);
        light;
        
        waitfor(vizRate);
    end

j = j+1;
end

function [t,phi,theta,psi] = controller_att(t, initial_x, initial_dxdt,theta_desired,theta_desired_dot)

    [t,phi]= ode45( @solv_att, t, [initial_x(1) initial_dxdt(1)] );
    [t,theta]=ode45( @solv_att, t, [initial_x(2) initial_dxdt(2)] );
    [t,psi]=ode45( @solv_att, t, [initial_x(3) initial_dxdt(3)] );
    
    phi(:,1) = phi(:,1) + theta_desired(1);
    phi(:,2) = phi(:,2) + theta_desired_dot(1);
    
    theta(:,1) = theta(:,1) + theta_desired(2);
    theta(:,2) = theta(:,2) + theta_desired_dot(2);
    
    psi(:,1) = psi(:,1) + theta_desired(3);
    psi(:,2) = psi(:,2) + theta_desired_dot(3);
end
 
function moments = dyn_inv_att(t, phi, theta, psi, theta_desired, theta_desired_dot, J)

    theta = [phi(:,1) + theta_desired(1),theta(:,1) + theta_desired(2),psi(:,1) + theta_desired(3)];
    thetadot = [phi(:,2) + theta_desired_dot(1),theta(:,2) + theta_desired_dot(2),psi(:,2) + theta_desired_dot(3)];
    thetadotdot = diff(thetadot)/0.001;

    thetadotdot( end+1, : ) = 0; 
    i=1;

    while (i<length(t)+1)
        th=theta(i,:);
        th_dot=thetadot(i,:);
        th_dotdot=thetadotdot(i,:);
        t1=[1, 0, -sin(th(2));
            0, cos(th(1)), sin(th(1))*cos(th(2));
            0, -sin(th(1)), cos(th(1))*sin(th(2));];
        t2=[0, 0, cos(th(2))*th_dot(2);
            0, sin(th(1))*th_dot(1), -sin(th(1))*sin(th(2))*th_dot(2)-cos(th(1))*cos(th(2))*th_dot(1);  
            0, cos(th(1))*th_dot(1), -sin(th(2))*cos(th(1))*th_dot(2)+sin(th(1))*cos(th(2))*th_dot(1);];

        omegadot=t1*transpose(th_dotdot)+t2*transpose(th_dot);
        omega=t1*transpose(th_dot);
        omega_hat=[0, -omega(3), omega(2);
                    omega(3), 0 , -omega(1);
                    -omega(2),omega(1), 0;];

        moment=omega_hat*J*omega+J*omegadot;
        moments(i,:)=transpose(moment);

        i=i+1; 
    end
end

function [t,x,y,z] = controller_pos(t, initial_x, initial_dxdt)

    [t,x]= ode45( @solv_pos, t, [initial_x(1) initial_dxdt(1)] );
    [t,y]=ode45( @solv_pos, t, [initial_x(2) initial_dxdt(2)] );
    [t,z]=ode45( @solv_pos, t, [initial_x(3) initial_dxdt(3)] );
end

function [Thrust, phi_desired,theta_desired,x_dot,y_dot] = dyn_inv_pos( x_error,y_error,z_error, desired_pos, desired_pos_dot) 

    pos = [x_error(:,1) + desired_pos(1),y_error(:,1) + desired_pos(2),z_error(:,1) + desired_pos(3)];
    pos_dot = [x_error(:,2) + desired_pos_dot(1),y_error(:,2) + desired_pos_dot(2),z_error(:,2) + desired_pos_dot(3)];
    pos_dotdot = diff(pos_dot)/0.5;

    i=1;
    m = 2;
    g = 9.81;
    
        
    Thrust = m*(sqrt((pos_dotdot(end,1))^2 + (pos_dotdot(end,2))^2 + (g - pos_dotdot(end,3))^2)) ;
    x_dot = -m*pos_dotdot(end,1)/Thrust ;
    y_dot = -m*pos_dotdot(end,2)/Thrust ;
    
    phi_desired = asin(-(y_dot));
    theta_desired = asin(x_dot/cos(phi_desired)) ;

    i=i+1; 
    
end


function dxdt=solv_pos(t,x)
    wn = 8.5;
    zeta = 0.5;
    dxdt_1 = x(2);
    dxdt_2 = -2*zeta*wn*x(2) - wn^2*x(1);
    dxdt=[dxdt_1; dxdt_2];
end


function dxdt=solv_att(t,x)
    wn = 10;
    zeta = 0.8;
    dxdt_1 = x(2);
    dxdt_2 = -2*zeta*wn*x(2) - wn^2*x(1);
    dxdt=[dxdt_1; dxdt_2];
end
 













