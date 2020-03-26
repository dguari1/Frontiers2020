T = 1/1000;
delay = 40/1000;

N = 30000;

INPUT_POS=load('C:\Users\GUARIND\Dropbox (Personal)\PhD\Paper_JointStiff_Movement\Code\position_input.mat');
input_position=INPUT_POS.position_input;
velocity=ddt_real(input_position)*(1/T);
acceleration=ddt_real(velocity)*(1/T);

K = 40;
B = 1.5;
I = 0.02;

Intrinsic_torque = K*input_position(1:N+100) + B*velocity(1:N+100) + I*acceleration(1:N+100);
Intrinsic_torque = Intrinsic_torque(101:N+100);


tau=ceil(delay/T); %discrete delay
del_velocity=[zeros(tau,1);velocity(1:(N+100)-tau)];  %delayed velocity
nl_del_velocity = del_velocity;
nl_del_velocity(nl_del_velocity<0)=0;

gainn= 0:0.5:20;

w = 16.1;
z = 0.65;

max_iteration = 10;
id_tolerance = 0.01;
tqI = zeros(N,length(gainn));
tqR = zeros(N,length(gainn));

mReflex=nlbl;
set(mReflex,'idMethod','sls', 'displayFlag',false,'nIterMax',10,'nLagLE',30, 'maxOrderNLE',6);


for k =1 : length(gainn)
    G = -gainn(k);
    Reflex_torque = lsim(tf(G*w*w,[1,2*z*w,w*w]), nl_del_velocity, (0:N-1+100)*T);

    
    Reflex_torque = Reflex_torque(101:N+100);


    RtoI(k,1) = var(Reflex_torque)/var(Intrinsic_torque);

    Total_torque = Intrinsic_torque+Reflex_torque;

    
    rip=0.05; [b,a] = cheby1(8,rip,0.8/10);Noise = filtfilt(b,a,randn(N,1)); %Gaussian noise, 50Hz bandwidth
    power_noise = sum(Noise.^2);
    power_signal = sum((Total_torque).^2);
    snr=15;
    Total_noise = Noise*sqrt((power_signal/(10^(snr/10)))/power_noise);

    Torque = Total_torque + Total_noise;

    
    tqI_temp = zeros(N/10,1);
    tqR_temp = zeros(N/10,1);
    
    VAFbest = -1000;
    for m = 1 : max_iteration
        tqI_temp = decimate(Torque,10) - tqR_temp;
        
        [Params_int,Ps_int,tqI_temp] = my_rivbjmiso(tqI_temp,[decimate(input_position(101:N+100),10) decimate(velocity(101:N+100),10) decimate(acceleration(101:N+100),10)],0,ones(1,3));

        tqR_temp = decimate(Torque,10)-tqI_temp;
        zv=cat(2,nldat(decimate(velocity(101:N+100),10),'domainIncr',1/100),nldat(tqR_temp,'domainIncr',1/100));
        mReflex=nlident(mReflex, zv);
        tqR_temp= double(nlsim(mReflex,nldat(decimate(velocity(101:N+100),10),'domainIncr',1/100)));
        
        tqT_temp = tqI_temp + tqR_temp;
        
        VAFT = VAFnl(decimate(Torque,10),tqT_temp);
        
        if (VAFT-VAFbest) < id_tolerance
            break
        else
            VAFbest=VAFT;
        end
    end
    
    
    [Params_int,Ps_int,est] = my_rivbjmiso(decimate(Torque,10) - tqR_temp,[decimate(input_position(101:N+100),10) decimate(velocity(101:N+100),10) decimate(acceleration(101:N+100),10)],0,ones(1,3));
    
    Elas_riv(k,1) = Params_int(1);
    Visc_riv(k,1) = Params_int(2);
    Iner_riv(k,1) = Params_int(3);
    
    Elas_var_riv(k,1) = sqrt(Ps_int(1,1));
    Visc_var_riv(k,1) = sqrt(Ps_int(2,2));
    Iner_var_riv(k,1) = sqrt(Ps_int(3,3));
    
    
    VAF_total(k,1) = VAFnl(decimate(Total_torque,10),est);
    VAF_int(k,1) = VAFnl(decimate(Intrinsic_torque,10),est);
   
end
%%
figure;
subplot(1,3,1)
errorbar(RtoI, Elas_riv, 2*Elas_var_riv)
hold on
plot(RtoI, K*ones(length(RtoI),1),'Color','red') 
subplot(1,3,2)
errorbar(RtoI, Visc_riv, 2*Visc_var_riv)
hold on
plot(RtoI, B*ones(length(RtoI),1),'Color','red') 
subplot(1,3,3)
errorbar(RtoI, Iner_riv, 2*Iner_var_riv)
hold on
plot(RtoI, I*ones(length(RtoI),1),'Color','red') 