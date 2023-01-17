close all

%% Initialisation----------------------------------------------------------
SDF2019aRead

% Time vector
t=data(:,1)-data(1,1);
dt=t(2);

% Standard deviation (see Tables)
sigma_x = 0.5; % m
sigma_y = 0.5; % m
% Acceleration (see Tables)
a = 0.005*1e-3; % m/s^2
% Position and velocity (due to unmodelled acceleration effects)
ds = 0.5*a*dt^2;
dv = a*dt;

u0=0;
v0=0;

%% MODEL INITIALIZATION---------------------------------------------------
% M: system transition matrix
M = [1,0,dt,1;0,1,0,dt;0,0,1,0;0,0,0,1];

% K: Measurement matrix
K = [1,0,0,0;0,1,0,0];

% S_xi: Process noise covariance matrix
S_xi = [ds^2,0,0,0;0,ds^2,0,0;0,0,dv^2,0;0,0,0,dv^2];

% S_eps: Measurement noise covariance matrix
S_eps= [sigma_x^2,0;0,sigma_y^2];

%% FILTER INITIALIZATION---------------------------------------------------
% x: State vector allocation - use the command zeros to initialize a
% matrix (time along columns)
x = zeros(4,length(data));
% x0: initial state vector
x(:,1) = [data(1,2),data(1,3),0,0];

% S: Covariance matrix allocation - use the command zeros to initialize a
% matrix (time along 3rd dimension)
S = zeros(size(x,1),size(x,1),length(data));

% S0: State Covariance matrix (initial value)
S (:,:,1)= [sigma_x^2,0,0,0;0,sigma_y^2,0,0;0,0,dv^2,0;0,0,0,dv^2];

% Measurements - definition of y matrix (remove time):
y = [data(:,2)';data(:,3)'];

%% Linear KF algorithm-----------------------------------------------------
for k = 2:length(data)
    %% Time Update ========================================================
    % Equation 1: State propagation
    xa =  M*x(:,k-1);
    % Equation 2: State covariance propagation
    Sa =  M*S(:,:,k-1)*M'+S_xi;
    Yzz=(K*Sa*K'+S_eps)^(-1);
    %% Measurement Update==================================================
    % Equation 3: Kalman Gain
    G = Sa*K'*Yzz ;
    % Equation 4: Innovation
    nu= y(:,k)-K*xa ;
    % Equation 5: State update
    x(:,k) = xa+G*nu ;
    % Equation 6: Error covariance update
    S(:,:,k) = Sa-G*K*Sa;
    
    %% Compute Chi-squared statistics======================================
    %Information matrix of innovations - notice you are computing Yzz
    %twice. You can compute Yzz before Equation 3 to make the code more
    %efficient.
    %Yzz=(K*Sa*K'+S_eps)^(-1);
    ChiStat(k)=nu'*Yzz*nu;
end



%% Plots-------------------------------------------------------------------
%Plots: time-histories-----------------------------------------------------
figure
subplot(2,1,1)
plot(t,data(:,2),t,x(1,:))
legend('x Measured','x Filtered')
title('Measured and Filtered measures vs time')
xlabel('Time [s]')
ylabel('Measured x [m]')

subplot(2,1,2)
plot(t,data(:,3),t,x(2,:))
legend('y Measured','y Filtered')
title('Measured and Filtered measures vs time')
xlabel('Time [s]')
ylabel('Measured x [m]')


% Plots: Rover LOCUS-----------------------------------------------------
figure
plot(y(1,:),y(2,:),x(1,:),x(2,:))
legend('Measured','Filtered')
title('Rover locus [X,Y]')
xlabel('X coordinate [m]')
ylabel('Y coordinate [m]')


% Plots: Chi-Squared statistics--------------------------------------------
figure;
plot(t,ChiStat)
xlabel('Time [s]')
ylabel('\chi^2')
title('\chi^2 statistics about innovations')
