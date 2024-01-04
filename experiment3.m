rng(2023);
Q_t1 =  5e6;
Q_t2 = 4e6;
sigma1 = 15;
sigma2 = 60;
t1 = 80;
t2 = 160;
t = 0:300;

 V1 = (t+sigma2-t1>0)*Q_t1.*(t+sigma2-t1)./(sigma2^2).*exp(-(t+sigma2-t1).^2/(2*sigma2^2));
 V2 = (t+sigma2-t2>0)*Q_t2.*(t+sigma2-t2)./(sigma2^2).*exp(-(t+sigma2-t2).^2/(2*sigma2^2));


lam = 1.1*max(max(V1),max(V2));
V = V1+V2;
plot(t,V1/lam);
hold on;
plot(t,V2/lam);
title('Waveform before superposition occurs');
xlabel('Time（ns）');
ylabel('Amplified signal(V)');
legend('Q_t_1','Q_t_2');
axis([0,300,0,1.2]);
figure;
plot(t,V/lam);
title('Waveform after superposition occurs');
xlabel('Time（ns）');
ylabel('Amplified signal(V)');
axis([0,300,0,1.2]);