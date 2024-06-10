clear all
%%% Path of Spectral function
addpath('C:/Users/fukas/Dropbox/git/spectral')

spectral_spec=1;
common_alpha_spec=0;
alpha0_param=1;
lambda_param=1e-7;
D=4;


Method=-1;
run Main_function.m
c0_temp0=c0;
iter_info00=iter_info;

if 1==1
Method=-2;
run Main_function.m
c0_temp00=c0;
iter_info000=iter_info;

Method=0;
run Main_function.m
c0_temp00=c0;
iter_info0=iter_info;

Method=1;
run Main_function.m
c0_temp1=c0;
iter_info1=iter_info;

Method=2;
run Main_function.m
c0_temp2=c0;
iter_info2=iter_info;
end
