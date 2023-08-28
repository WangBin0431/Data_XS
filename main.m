rng(2023);
Position = struct('retain',[],'remove',[]);
Position = repmat(Position,[1,1000]);
DATA = struct('s',0,'T',0,'sigma_p',0,'x0',0,'y0',0,'Position',Position );
DATA = repmat(DATA,[1,80]);
DATA_num = 1;


rate = 8;   

s = 1; 
K = 10;
D_p = 5;   
sigma_p = D_p*rate/6;   
var_value = 10:10:200;  
var_name = 'K';     
var_num = size(var_value,2);    
expe_data;

s = 1;    
K = 10;   
D_p = 5;     
sigma_p = D_p*rate/6;   
var_value = 0.5:0.5:10;
var_name = 's';     
var_num = size(var_value,2);    
expe_data;

s = 1;    
K = 100;   
D_p = 5;     
sigma_p = D_p*rate/6;   
var_value = 0.5:0.5:10;  
var_name = 's';     
var_num = size(var_value,2);    
expe_data;

s = 1; 
K = 100;
D_p = 5;   
sigma_p = D_p*rate/6;   
var_value = 2:0.5:11.5;  
var_name = 'p';     
var_num = size(var_value,2);    
expe_data;



save_name = 'DATA.mat';
save(save_name,'DATA');