theta_desired = [0;0;0];
theta = deg2rad(2 * deviation * rand(3, 1) - deviation);
e = theta - theta_desired;

theta_desired_dot = [0;0;0];
theta_dot = [0;0;0];
e_dot = theta_dot - theta_desired_dot;

t = 0:0.01:10; 

J = [0.01, 0, 0; 
     0, 0.05, 0; 
     0, 0, 0.001];

initial_x = e;
initial_dxdt = e_dot;

        [t,phi,theta,psi] = controller_att( t, initial_x, initial_dxdt) ;

        moments = dyn_inv_att(t, phi, theta, psi, theta_desired, theta_desired_dot,J) ;
        
        realTheta_1 = phi(:,2);
        realTheta_2 = theta(:,2);
        realTheta_3 = psi(:,2);
                   
t = tiledlayout(1,3); % Requires R2019b or later
nexttile
plot(realTheta_1)
nexttile
plot(realTheta_2)
nexttile
plot(realTheta_3)

t.Padding = 'compact';
t.TileSpacing = 'compact';

    
    
    
    function [t,phi,theta,psi] = controller_att(t, initial_x, initial_dxdt)

    [t,phi]= ode45( @solv_att, t, [initial_x(1) initial_dxdt(1)] );
    [t,theta]=ode45( @solv_att, t, [initial_x(2) initial_dxdt(2)] );
    [t,psi]=ode45( @solv_att, t, [initial_x(3) initial_dxdt(3)] );
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

    function dxdt=solv_att(t,x)
    wn = 10;
    zeta = 0.2;
    dxdt_1 = x(2);
    dxdt_2 = -2*zeta*wn*x(2) - wn^2*x(1);
    dxdt=[dxdt_1; dxdt_2];
end