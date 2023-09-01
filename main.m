% The number of photons emitted from the detected target to the detector in 100ns is s.
s = 1; 
%Total detection time in 100ns
T = 10;   
%Photon diffusion diameter, in number of anode spacers
D_r = 5;  

var_value = 10:10:200;  
var_name = 'T';   
experiment(s, T, D_r, var_value, var_name);
var_value = 0.5:0.5:10;  
var_name = 's';   
experiment(s, T, D_r, var_value, var_name);

T = 100;   
var_value = 0.5:0.5:10;  
var_name = 's';   
experiment(s, T, D_r, var_value, var_name);
var_value = 2:0.5:11.5;  
var_name = 'D';   
experiment(s, T, D_r, var_value, var_name);