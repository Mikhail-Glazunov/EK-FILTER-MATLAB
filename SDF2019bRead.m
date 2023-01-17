%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aerospace Navigation and Sensor assignment
% AY 2019-20
% Cranfield University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

%Import SDF2019b.txt file
filename='SDF2019b.txt';
data=dlmread(filename,' ',3,0);


%Plot data
h=figure
plot(data(:,2),data(:,3))
title('Measured coordinates of the rover')
xlabel('X coordinate [m]')
ylabel('Y coordinate [m]')

%Plot data vs time
h=figure
plot(data(:,1),data(:,2:3))
title('Measured coordinates vs time')
xlabel('Time [s]')
ylabel('Measured position [m]')

%Plot data vs time
h=figure
subplot(2,1,1),plot(data(:,1),data(:,4))
title('Measured range vs time')
xlabel('Time [s]')
ylabel('Measured range [m]')

subplot(2,1,2),plot(data(:,1),data(:,5))
title('Measured angle vs time')
xlabel('Time [s]')
ylabel('Angle [rad]')

