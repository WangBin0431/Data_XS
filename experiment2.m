rng(2023);
t1 = 50;
dt = 170;
t2 = t1+dt;

sigma_Q = 80;
rate = 0.2:0.05:5;
Q_total_1= 6e6;
d=1;
x0 = 20.3856346;
sigma_e=5/6;
num =1000;

t = 0:300;
po = zeros(1,size(rate,2));
for j = 1:size(rate,2)
    Q_total_2= Q_total_1*rate(j);

    x1 = zeros(1,num);
    x2 = zeros(1,num);

    x_position = zeros(1,num);
    for ex =1:num
        x1(ex) = normrnd(x0,sigma_e);
        x2(ex) = normrnd(x0,sigma_e);
        Qx = zeros(1,40);
        Qy = zeros(1,40);

        Fuzhi_x1 = zeros(1,40);
        Fuzhi_x2 = zeros(1,40);

        for i3 = 1:40
            Fuzhi_x1(i3) = Q_total_1/((2*pi)^(1/2)*sigma_e)*exp(-(i3-x1(ex))^2/(2*sigma_e^2))*d;
            Fuzhi_x2(i3) = Q_total_2/((2*pi)^(1/2)*sigma_e)*exp(-(i3-x2(ex))^2/(2*sigma_e^2))*d;

        end

        for i3 = 10:30
            wave_x = (t+sigma_Q-t1>=0)*Fuzhi_x1(i3).*(t+sigma_Q-t1)./(sigma_Q^2).*exp(-(t+sigma_Q-t1).^2/(2*sigma_Q^2))+...
                (t+sigma_Q-t2>=0)*Fuzhi_x2(i3).*(t+sigma_Q-t2)./(sigma_Q^2).*exp(-(t+sigma_Q-t2).^2/(2*sigma_Q^2));

            Qx(i3) = max(wave_x);
        end


        max_Qx = find(Qx == max(Qx));


        x_position(ex) = Qx(max_Qx-2:max_Qx+2) * (max_Qx-2:max_Qx+2)' / sum(Qx(max_Qx-2:max_Qx+2));

    end


    
    
    x1 = abs(x1-x0);
    x2 = abs(x2-x0);
    x3 = abs(x_position-x0);
    xsum = [x1,x2]; 
    
    x1_ = sum(x1)/size(x1,2);
    x2_ = sum(x2)/size(x2,2);
    x3_ = sum(x3)/size(x3,2);
    xsum_ = sum(xsum)/size(xsum,2);
    
    p11(j) = ((x1-x1_)*(x1-x1_)'/(size(x1,2)-1))^(1/2);
    p21(j) = ((x2-x2_)*(x2-x2_)'/(size(x2,2)-1))^(1/2);
    p31(j) = ((x3-x3_)*(x3-x3_)'/(size(x3,2)-1))^(1/2);
    psum1(j) = ((xsum-xsum_)*(xsum-xsum_)'/(size(xsum,2)-1))^(1/2);
    
    
    p1(j) = x1_;
    p2(j) = x2_;
    p3(j) = x3_;
    psum(j) = xsum_;

end


plot(rate,p3);
hold on;
plot(rate,psum);
title(['Interval time ',num2str(t2-t1),'ns']);
xlabel('Q_t_o_t_a_l(2)/Q_t_o_t_a_l(1)');
ylabel('Mean value of error (Anode width)');
legend('Aliased data','Separate data');
axis([0,size(rate,2)/20,0.8*min(min(p3),min(psum)),1.2*max(max(p3),max(psum))]);

figure;

plot(rate,p31);
hold on;
plot(rate,psum1);  
title(['Interval time ',num2str(t2-t1),'ns']);
xlabel('Q_t_o_t_a_l(2)/Q_t_o_t_a_l(1)');
ylabel('Standard deviation of error (Anode width)');
legend('Aliased data','Separate data');
axis([0,size(rate,2)/20,0.8*min(min(p31),min(psum1)),1.2*max(max(p31),max(psum1))]);