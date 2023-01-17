clc 
clear all
close all

%Kalmann Filter Implementation with only ONE landmark:
 i=3;
%% data loading
load('my_input.mat');
% v input linear velocity (m/s)
% om input angular velocity (rad/s)
% t time vector (s)

load('my_measurements.mat');
% r range measurment (m)
% b bearing measurments (rad)
% l landmark position

T=t(2)-t(1); %sample period

r=double(r);
b=double(b);
l=double(l);



%% Matrix Definition

% f_M: transition function (non-linear)
f_M = @(x,v,om) [x(1,:) + T*cos(x(3,:)).*v ; 
                 x(2,:) + T*sin(x(3,:)).*v ;
                 x(3,:) + T*om];

%S_xi: Process noise covariance matrix 
    
    sigma_v=0.008^0.5; %
    sigma_om=10^0.5; %
    
    S_xi=diag([sigma_v^2,sigma_om^2]);

% f_Ki: Measurement function for the landmark i 
    f_Ki = @(x,i) [sqrt( (l(i,1)-x(1,:) ).^2 + ( l(i,2)- x(2,:) ).^2 );
                   wrapToPi( atan2( l(i,2)- x(2,:) , l(i,1)- x(1,:) )) - x(3,:) ];

% Measurments error matrix (2 measurments x 2 measurments)      
    sigma_r=0.01^0.5; 
    sigma_b=0.01^0.5;
    S_eps=diag([sigma_r^2 sigma_b^2]);


%% algorithm

%initialisation of the state vector (given)
    
    x=zeros(3,length(b)); 
    
    x0=50; 
    y0=0; 
    theta0=pi/2; 
    x(:,1)=[x0 y0 theta0]';

%State Covariance matrix (initial value)
    S = zeros(3,3,length(b));
    
    dx=std(v)*T; % how to estimate it ?
    dy=std(v)*T; % how to estimate it ?
    dtheta=std(om)*T; % how to estimate it ?
    
    S(:,:,1) = [dx 0 0;
                0 dy 0; 
                0 0  dtheta]; %how to estimate it?

% Measurements initialisation
    y=[r(:,i)' ; b(:,i)'];

for k=2:501
    
    %% Time Update ========================================================
    % Equation 1: State propagation (non-linear equation)
    
    xa = f_M( x(:,k-1),v(k-1),om(k-1) )
    
    % Equation 2: State covariance propagation (linearised equation)
       
        %Evaluation of the M-matrix at k-1
        M=[1  0 -T*v(k-1)*sin(x(3,k-1));
           0  1  T*v(k-1)*cos(x(3,k-1));
           0  0  1];

       %Evaluation of the Gamma-matrix at k-1
        Gamma=[T*cos(x(3,k-1))  0  ;
               T*sin(x(3,k-1))  0  ;
               0                T];
          

        Sa =  M*S(:,:,k-1)*M'  + Gamma*S_xi*Gamma'; %Q=S_xi ??

    
    %% Measurement Update==================================================
    
    %Equation 3:prediction of the measurment (non-linear equation)
    ya = f_Ki( x(:,k-1),i ) ; 
    
    %Equation 3': Jacobian matrix (linearisation of the measurement model)
    
        %Construction of the K-matrix
        %K=[]      
        %Construction of the K_i matrices (2 x 3)
        A=( l(i,1)-x(1,k) ).^2 + ( l(i,2)- x(2,k) ).^2; 
        K11= (x(1,k)-l(i,1))./sqrt(A);
        K12= (x(2,k)-l(i,2))./sqrt(A);
        K13= 0;
        K21=-(x(2,k)-l(i,2))./A;
        K22= (x(1,k)-l(i,1))./A;
        K23= -1;

        K=[K11 K12 K13;
             K21 K22 K23];
        %Assembly of the K-matrix (16 x 3)
        %K=[K;K_i];
          
    % Equation 4: Kalman Gain
    G = Sa*K'*(K*Sa*K'+S_eps)^(-1);
    % Equation 5: Innovation
    nu= y(:,k) - f_Ki(xa,i) ;
    % Equation 6: State update
    x(:,k) = xa + G*nu;
    % Equation 7: Error covariance update
    S(:,:,k) = Sa - (G * K * Sa);
    

    %Chi statistique
    Yzz=(K*Sa*K'+S_eps)^(-1);
    ChiStat(k)=nu'*Yzz*nu;

end


%% Results Plotting
y_post=f_Ki(x(:,:),i);

figure(1);
plot(t,y(1,:),t,y_post(1,:),'.');
xlabel('t [s]');
ylabel('range r [m]');
legend('Measured','Filtered');
title('r');

figure(2);
plot(t,y(2,:),t,y_post(2,:),'.');
xlabel('t [s]');
ylabel('angle b [rad]');
legend('Measured','Filtered');
title('b');

figure(3);
plot(t,v);
xlabel('t [s]');
ylabel('v (m/s)');
title('v input');

figure(4);
plot(t,om);
xlabel('t [s]');
ylabel('\omega (rad/s)');
title('\omega input');

figure(5);
plot(t(1:length(ChiStat)),ChiStat),
xlabel('Time [s]')
ylabel('\chi^2')
title('\chi^2 statistics about innovations')
