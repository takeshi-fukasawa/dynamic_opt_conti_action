clear
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
