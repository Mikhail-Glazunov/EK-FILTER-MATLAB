close all
%% Initialisation-------------------------------------------------------
SDF2019bRead

% Time vector
t=data(:,1)-data(1,1);
dt=t(2);

% Standard deviation from datasheet
sigma_x = .5; % m/s
sigma_y = .5; % m/s
sigma_r =  1;
sigma_theta= 1e-4;
% Acceleration from datasheet
a = 1E-6; % m/s^2
%Position and velocity
ds = 0.5*a*dt^2;
dv = a*dt;


u0=0;
v0=0;


%% MODEL INITIALIZATION---------------------------------------------------
% M: system transition matrix
M = [1 0 dt 0;
    0 1 0 dt;
    0 0 1 0;
    0 0 0 1];

% f_K: Measurement matrix (non-linear) - you can use an anonymous function
% with argument input @(x)
% y - f_K(x)
% y - [x y r theta]
% y - f_R([x y xdot ydot])
f_K = @(x) [x(1);x(2);sqrt((x(1)-u0).^2+(x(2)-v0).^2);wrapToPi(atan2((x(2)-v0),(x(1)-u0)))];

% S_xi: Process noise covariance matrix
S_xi = [ds^2 0 0 0;
     0 ds^2 0 0;
     0 0 dv^2 0;
     0 0 0 dv^2];

% S_eps: Measurement noise covariance matrix
S_eps= [sigma_x^2,0,0,0;0,sigma_y^2,0,0;0,0,sigma_r^2,0;0,0,0,sigma_theta^2];

%% FILTER INITIALIZATION---------------------------------------------------
% x: State vector allocation - use the command zeros to initialize a
% matrix (time along columns)
x = zeros(4,length(data));
%initial state vector
x(:,1) = [data(1,2) data(1,3) 0 0];

% S: Covariance matrix allocation - use the command zeros to initialize a
% matrix (time along 3rd dimension)
S = zeros(size(x,1),size(x,1),length(data));

% S0: State Covariance matrix (initial value)
S(:,:,1) = [sigma_x^2 0 0 0;
    0 sigma_y^2 0 0;
    0 0 dv^2 0;
    0 0 0 dv^2];

% Measurements initialisation
y = [data(:,2)';data(:,3)';data(:,4)';data(:,5)'];

%% EKF Filter algorithm-----------------------------------------------------
for k = 2:length(data)
    %% Time Update ========================================================
    % Equation 1: State propagation
    xa = M*x(:,k-1);
    % Equation 2: State covariance propagation
    Sa = (M * S(:,:,k-1) * (M.')) + S_xi;

    
    %% Measurement Update==================================================
    % Compute the non-linear prediction of the measurement (useful to
    % compute the next step - Jacobian):
    ya= f_K(xa);
    
    
    % Equation 3: Jacobian matrix (linearization of the measurement model)
    % (m x n) -> 4 measurements x 4 states 
    K = [1,0,0,0;0,1,0,0;cos(ya(4)),sin(ya(4)),0,0;-sin(ya(4))/ya(3),cos(ya(4))/ya(3),0,0];
    
    % Equation 4: Kalman Gain
    G = Sa*K'*(K*Sa*K'+S_eps)^(-1);
    % Equation 5: Innovation
    nu=  y(:,k)-ya;
    % Equation 6: State update
    x(:,k) = xa + G*nu;
    % Equation 7: Error covariance update
    S(:,:,k) = Sa - (G * K * Sa);
    
    %% Compute Chi-squared statistics======================================
    %Information matrix of innovations - notice you are computing Yzz
    %twice. You can compute Yzz before Equation 3 to make the code more
    %efficient.
    Yzz=(K*Sa*K'+S_eps)^(-1);
    ChiStat(k)=nu'*Yzz*nu;
end

%% POST-PROCESSING
%Compute the filtered measurements (using filtered states and the full non-linear
%measurement model):
y_post= f_K(x);

%Plots: time-histories-----------------------------------------------------
figure
subplot(2,2,1)
plot(t,data(:,2),t,y_post(1,:))
legend('x Measured','x Filtered')
title('Measured and Filtered measures vs time')
xlabel('Time [s]')
ylabel('Measured x [m]')

subplot(2,2,2)
plot(t,data(:,3),t,y_post(2,:))
legend('y Measured','y Filtered')
title('Measured and Filtered measures vs time')
xlabel('Time [s]')
ylabel('Measured x [m]')

subplot(2,2,3)
plot(t,data(:,4),t,y_post(3,:))
legend('r Measured','r Filtered')
title('Measured and Filtered range vs time')
xlabel('Time [s]')
ylabel('Measured range r [m]')

subplot(2,2,4)
plot(t,data(:,5),t,y_post(4,:))
title('Measured and Filtered bearing angle vs time')
legend('\theta Measured','\theta Filtered')
xlabel('Time [s]')
ylabel('Measured \theta [rad]')


% Plots: Rover LOCUS-----------------------------------------------------
figure
hold on 
plot(x(1,:),x(2,:))
plot(Xa(1,:),Xa(2,:))
legend('Filtered','Measured')
title('Rover locus [X,Y]')
xlabel('X coordinate [m]')
ylabel('Y coordinate [m]')
hold off

% Plots: Chi-Squared statistics--------------------------------------------
figure;
plot(t,ChiStat)
xlabel('Time [s]')
ylabel('\chi^2')
title('\chi^2 statistics about innovations')