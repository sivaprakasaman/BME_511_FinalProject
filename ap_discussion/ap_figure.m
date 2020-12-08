%show alternating polarity/env/tfs relationship
close all

fs = 1000;
time = 2.75;
start = .75;
t = start:1/fs:time;
cr = 15;
mod = 1;

p = (sin(2*pi*mod*t)+1).*(sin(2*pi*cr*t)); 
env = sin(2*pi*mod*t)+1;
hold on
plot(t,p,'LineWidth',1.5)
plot(t,-p,'LineWidth',1.5)
plot(t,env,'k','LineWidth',1.5);
hold off
legend('Positive','Negative','Envelope')
xlabel('Time (s)')
ylabel('Amplitude');