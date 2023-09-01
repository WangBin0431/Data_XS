function [] = experiment(s, T, D_r, var_value, var_name) 

rng(2023);
% Quantum efficiency
n = 0.2;   
%Position of the object to be measured
xa = 100+rand()*100;      
ya = 100+rand()*100; 
%Number of pixels in x-direction
X_max = 320;   
%Number of pixels in y-direction
Y_max = 320;    
%Each anode interval is equivalent to 8 pixels
rate = 8;  
%Electron cloud diameter in number of anode spacers
D_e = 5; 
%Dead time in 100ns
t_d = 2;    
%Standard deviation of the photon dispersion spot. Unit: pixel
sigma_r = D_r*rate/6;   
%Standard deviation of the electron cloud, in units of the number of anode spacings
sigma_e = D_e/6;    
%Width of anode
d = 0.1;    
P = 1;  
% Variance of Gaussian waveforms
sigma_Q = 0.8;  

var_num = size(var_value,2);    
r_u = zeros([1,var_num]);  
r_s = zeros([1,var_num]);   
r_u2 = zeros([1,var_num]);  
r_s2 = zeros([1,var_num]);  
exp_num = 1000;  
poi_num = 8;    

for var_i = 1:var_num
    if var_name == 'T'
        T = var_value(var_i);
    elseif var_name == 'D'
        sigma_r=var_value(var_i)*rate/6;
    elseif var_name == 's'
        s = var_value(var_i);
    end          
    r_1 = [];
    r_2 = [];
    for exp_i =1:exp_num        
        Photon_time = [];
        Time = 0;
        %Get the arrival time of the photons
        while Time<T
            Arr_time = Time + exprnd(1/(n*s));
            if Arr_time>T
                break;
            end
            Photon_time = [Photon_time,Arr_time];
            Time = Arr_time;
        end
        Photon_num = size(Photon_time,2);
        if Photon_num == 0
            continue ;
        end
        Photon_time_delta = Photon_time - [Photon_time(1),Photon_time(1:Photon_num-1)];
        %Get the arrival position of the photon
        Photon_x = normrnd(xa,sigma_r,[1,Photon_num]);
        Photon_y = normrnd(ya,sigma_r,[1,Photon_num]);
 
        no_mix_position = [];
        mix_position = [];
        ii0 = 1;
        while ii0 <= Photon_num    
            ii1 =ii0+1;
            %Determine the aliasing of each photon
            mix_num = 0;
            while 1
                if ii1 > Photon_num
                    break;
                elseif Photon_time_delta(ii1)>t_d
                    break;
                else
                    ii1 =ii1+1;
                    mix_num = mix_num + 1;
                end
            end

            Qx = zeros(1,40);
            Qy = zeros(1,40);
            if mix_num == 0
                %Processing of non-aliased data
                Q_total = (poissrnd(poi_num)+0.01)*2.5e6;
                x_r = P*(1:X_max/rate);
                y_r = P*(1:Y_max/rate);
                fx = Q_total/((2*pi)^(1/2)*sigma_e)*exp(-(x_r-Photon_x(ii0)/rate).^2/(2*sigma_e^2));
                fy = Q_total/((2*pi)^(1/2)*sigma_e)*exp(-(y_r-Photon_y(ii0)/rate).^2/(2*sigma_e^2));     
                Qx = fx*d;
                Qy = fy*d;
                max_Qx = find(Qx == max(Qx));
                max_Qy = find(Qy == max(Qy));
                max_Qx = max_Qx(1);
                max_Qy = max_Qy(1);
                x_position = round(rate * P * Qx(max_Qx-2:max_Qx+2) * (max_Qx-2:max_Qx+2)' / sum(Qx(max_Qx-2:max_Qx+2)));
                y_position = round(rate * P * Qy(max_Qy-2:max_Qy+2) * (max_Qy-2:max_Qy+2)' / sum(Qy(max_Qy-2:max_Qy+2)));
                no_mix_position = [no_mix_position;x_position,y_position];
 
            else
                %Processing of aliased data
                Q_total = (poissrnd(poi_num,1,mix_num+1)+0.01)*2.5e6;
                Q_max_num = find(Q_total == max(Q_total));
                f = 0 ;
                Q_t = Photon_time(ii0:ii0+mix_num);
                Peak_maxt = max(Q_t);
                Peak_mint = min(Q_t);
                t= Peak_mint:0.01:Peak_maxt;
                x_r = P*(1:X_max/rate);
                y_r = P*(1:Y_max/rate);
                for i3 = 1:X_max/rate                   
                    Peak_x = Q_total./((2*pi)^(1/2)*sigma_e).*exp(-(i3*P-Photon_x(ii0:ii0+mix_num)/rate).^2/(2*sigma_e^2))*d;
                    Peak_y = Q_total./((2*pi)^(1/2)*sigma_e).*exp(-(i3*P-Photon_y(ii0:ii0+mix_num)/rate).^2/(2*sigma_e^2))*d;
                    wave_x = 0;
                    wave_y = 0;
                    for i4 = 1:mix_num+1
                        wave_x = wave_x+Peak_x(i4)*(t+sigma_Q-Q_t(i4)>0).*...
                            (t+sigma_Q-Q_t(i4))./(sigma_Q^2).*exp(-(t+sigma_Q-Q_t(i4)).^2/(2*sigma_Q^2));
                        wave_y = wave_y+Peak_y(i4)*(t+sigma_Q-Q_t(i4)>0).*...
                            (t+sigma_Q-Q_t(i4))./(sigma_Q^2).*exp(-(t+sigma_Q-Q_t(i4)).^2/(2*sigma_Q^2));                       
                    end                    
                    Qx(i3) = max(wave_x);
                    Qy(i3) = max(wave_y);
                end
                
                max_Qx = find(Qx == max(Qx));
                max_Qy = find(Qy == max(Qy));
                x_position = round(rate * P * Qx(max_Qx-2:max_Qx+2) * (max_Qx-2:max_Qx+2)' / sum(Qx(max_Qx-2:max_Qx+2)));
                y_position = round(rate * P * Qy(max_Qy-2:max_Qy+2) * (max_Qy-2:max_Qy+2)' / sum(Qy(max_Qy-2:max_Qy+2)));
                mix_position = [mix_position;x_position,y_position];
            end
            ii0 = ii0 + mix_num + 1;
        end
        

        cal_2 = mean([no_mix_position;mix_position],1);
        cal_2_d = ((cal_2(1)-xa)^2+(cal_2(2)-ya)^2)^(1/2);
        r_2 = [r_2,cal_2_d];
        if size(no_mix_position,2)==0
            continue;
        end
        cal_1 = mean(no_mix_position,1);           
        cal_1_d = ((cal_1(1)-xa)^2+(cal_1(2)-ya)^2)^(1/2);
        r_1 = [r_1,cal_1_d];      
    end 

    r_u(var_i) = sum(r_1)/size(r_1,2);
    r_s(var_i) = ((r_1-r_u(var_i))*(r_1-r_u(var_i))'/(size(r_1,2)-1))^(1/2);
    r_u2(var_i) = sum(r_2)/size(r_2,2);
    r_s2(var_i) = ((r_2-r_u2(var_i))*(r_2-r_u2(var_i))'/(size(r_2,2)-1))^(1/2);

end

label_name = 'label';
if var_name == 'T'
    var_value = var_value/10;
    label_name = 'T(\mus)';
elseif var_name == 'D'
    label_name= 'D_r';
elseif var_name == 's'
    label_name = 's';
end
figure;
err_max = max(r_u)+max(r_s)+1;
errorbar(var_value,r_u,r_s,'red');
axis([var_value(1)-0.5,var_value(var_num)+1,0,1.2*err_max]);

hold on;
errorbar(var_value,r_u2,r_s2,'blue');
axis([var_value(1)-0.5,var_value(var_num)+1,0,1.2*err_max]);
xlabel(label_name);
ylabel('error(px)');
legend('removing','retaining');
end
