%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aerospace Navigation and Sensor Practical
% AY 2019-20
% Cranfield University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

%Import SDF2019a.txt file
filename='SDF2019a.txt';
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



