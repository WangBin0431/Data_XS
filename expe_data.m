
n = 0.2;    
x0 = 100+rand()*100;      
y0 = 100+rand()*100; 
X_max = 320;    
Y_max = 320;    
rate = 8;   


D_e = 5;    

t_d = 2;    


sigma_e = D_e/6;    
d = 0.1;    
P = 1;  

sigma_Q = 0.8;  



var_num = size(var_value,2);   
r_u = zeros([1,var_num]);  
r_s = zeros([1,var_num]);   
r_u2 = zeros([1,var_num]);  
r_s2 = zeros([1,var_num]);  
exp_num = 1000;  
poi_num = 8;    

for var_i = 1:var_num
    if var_name == 'K'
        K = var_value(var_i);
    elseif var_name == 'p'
        sigma_p=var_value(var_i)*rate/6;
    elseif var_name == 's'
        s = var_value(var_i);
    end
    DATA(DATA_num).s = s;
    DATA(DATA_num).T = 0.1*K;
    DATA(DATA_num).sigma_p = sigma_p;
    DATA(DATA_num).x0 = x0;
    DATA(DATA_num).y0 = y0;
    
    

    r_1 = [];
    r_2 = [];

    for exp_i =1:exp_num        
        Photon_time = [];
        Time = 0;

        while Time<K
            Arr_time = Time + exprnd(1/(n*s));
            if Arr_time>K
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

        Photon_x = normrnd(x0,sigma_p,[1,Photon_num]);
        Photon_y = normrnd(y0,sigma_p,[1,Photon_num]);

        no_mix_position = [];

        mix_position = [];
        ii0 = 1;

        while ii0 <= Photon_num    
            ii1 =ii0+1;

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
                Q_total = (poissrnd(poi_num,1,mix_num+1)+0.01)*2.5e6;
                Q_max_num = find(Q_total == max(Q_total));
                f = 0 ;
                Q_t = Photon_time(ii0:ii0+mix_num);
                Fuzhi_maxt = max(Q_t);
                Fuzhi_mint = min(Q_t);
                t= Fuzhi_mint:0.01:Fuzhi_maxt;

                Fuzhi_Qx_i = zeros(mix_num+1,X_max/rate);
                Fuzhi_Qy_i = zeros(mix_num+1,Y_max/rate);
                x_r = P*(1:X_max/rate);
                y_r = P*(1:Y_max/rate);
                for i3 = 1:X_max/rate
                    
                    Fuzhi_x = Q_total./((2*pi)^(1/2)*sigma_e).*exp(-(i3*P-Photon_x(ii0:ii0+mix_num)/rate).^2/(2*sigma_e^2))*d;
                    Fuzhi_y = Q_total./((2*pi)^(1/2)*sigma_e).*exp(-(i3*P-Photon_y(ii0:ii0+mix_num)/rate).^2/(2*sigma_e^2))*d;
                    wave_x = 0;
                    wave_y = 0;

                    for i4 = 1:mix_num+1
                        wave_x = wave_x+Fuzhi_x(i4)*(t+sigma_Q-Q_t(i4)>0).*...
                            (t+sigma_Q-Q_t(i4))./(sigma_Q^2).*exp(-(t+sigma_Q-Q_t(i4)).^2/(2*sigma_Q^2));
                        wave_y = wave_y+Fuzhi_y(i4)*(t+sigma_Q-Q_t(i4)>0).*...
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
        

        retain = [no_mix_position;mix_position];


        remove = no_mix_position;
        
        

        DATA(DATA_num).Position(exp_i).retain = retain;
        DATA(DATA_num).Position(exp_i).remove = remove;
        

        
    end 
    
    DATA_num = DATA_num +1;


end

