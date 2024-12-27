clearvars;close all;clc;
% contact matrix
M = [19.2 4.8 3.0 7.1 3.7 3.1 2.3 1.4 1.4;
     4.8 42.4 6.4 5.4 7.5 5.0 1.8 1.7 1.7;
     3.0 6.4 20.7 9.2 7.1 6.3 2.0 0.9 0.9;
     7.1 5.4 9.2 16.9 10.1 6.8 3.4 1.5 1.5;
     3.7 7.5 7.1 10.1 13.1 7.4 2.6 2.1 2.1;
     3.1 5.0 6.3 6.8 7.4 10.4 3.5 1.8 1.8;
     2.3 1.8 2.0 3.4 2.6 3.5 7.5 3.2 3.2;
     1.4 1.7 0.9 1.5 2.1 1.8 3.2 7.2 7.2;
     1.4 1.7 0.9 1.5 2.1 1.8 3.2 7.2 7.2];

load covidItalia28Oct.mat

figure(1)
subplot(2,3,1)
plot(timeVector,confirmed, 'linewidth', 2)      %totale casi confermati
xlabel('time');title('Total Confirmed cases')
set(gca,'FontSize',12)
subplot(2,3,2)