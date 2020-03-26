T = 1/1000;
delay = 40/1000;

N = 30000;

INPUT_POS=load('C:\Users\GUARIND\Dropbox (Personal)\PhD\Paper_JointStiff_Movement\Code\position_input.mat');
input_position=INPUT_POS.position_input;
velocity=ddt_real(input_position)*(1/T);
acceleration=ddt_real(velocity)*(1/T);


NOISE=load('C:\Users\GUARIND\Dropbox (Personal)\PhD\Paper_JointStiff_Movement\Code\Noise.mat');
out_noise=NOISE.Noise;


K = 40;
B = 0.5;
I = 0.02;

Intrinsic_torque = K*input_position(1:N+100) + B*velocity(1:N+100) + I*acceleration(1:N+100);
Intrinsic_torque = Intrinsic_torque(101:N+100);


tau=ceil(delay/T); %discrete delay
del_velocity=[zeros(tau,1);velocity(1:(N+100)-tau)];  %delayed velocity
nl_del_velocity = del_velocity;
nl_del_velocity(nl_del_velocity<0)=0;

gainn= 0:0.5:20;

G = -14.5;
w = 16.1;
z = 0.65;

for k =1 : length(gainn)
    G = -gainn(k);
    Reflex_torque = lsim(tf(G*w*w,[1,2*z*w,w*w]), nl_del_velocity, (0:N-1+100)*T);

    
    Reflex_torque = Reflex_torque(101:N+100);


    RtoI(k,1) = var(Reflex_torque)/var(Intrinsic_torque);

    Total_torque = Intrinsic_torque+Reflex_torque;


    %Noise=out_noise(1:N);
    %Noise = randn(N,1);
    rip=0.05; [b,a] = cheby1(8,rip,0.8/10);Noise = filtfilt(b,a,randn(N,1));
    power_noise = sum(Noise.^2);
    power_signal = sum((Total_torque).^2);
    snr=15;
    Total_noise = Noise*sqrt((power_signal/(10^(snr/10)))/power_noise);

    Torque = Total_torque + Total_noise;

    % figure;
    % subplot(3,1,1)
    % plot(Intrinsic_torque)
    % ylim([-20 20])
    % subplot(3,1,2)
    % plot(Reflex_torque)
    % ylim([-20 20])
    % subplot(3,1,3)
    % plot(Torque)
    % ylim([-20 20])

%     %estimation
%     pos = decimate(input_position(101:N+100),10);
%     tor = decimate(Torque,10);
% 
% 
%     nlag_intrinsic=ceil(.04/(1/100));
%     hIntrinsic = irf;
%     hIntrinsic= set(hIntrinsic,'nLags',nlag_intrinsic,'nSides',2);
% 
%     zp=cat(2,nldat(pos,'domainIncr',1/100),nldat(tor,'domainIncr',1/100));
%     hIntrinsic = nlident(hIntrinsic,zp);
%     %figure;plot((-4:4)*(1/100),hIntrinsic.dataSet)
% 
% 
%     tqI=nlsim(hIntrinsic,nldat(pos,'domainIncr',1/100));
%     Vafy(k,1)=VAFnl(decimate(Intrinsic_torque,10),tqI.dataSet);
    
%     [est,pars,noise_var,lk,COV]=TV_Bayes(decimate(Torque,10),[decimate(input_position(101:N+100),10) decimate(velocity(101:N+100),10) decimate(acceleration(101:N+100),10)]);
%     Elas(k,1) = pars(1);
%     Visc(k,1) = pars(2);
%     Iner(k,1) = pars(3);
%     
%     Elas_var(k,1) = sqrt(COV(1,1));
%     Visc_var(k,1) = sqrt(COV(2,2));
%     Iner_var(k,1) = sqrt(COV(3,3));
    
    [Params_int,Ps_int,est] = my_rivbjmiso(decimate(Torque,10),[decimate(input_position(101:N+100),10) decimate(velocity(101:N+100),10) decimate(acceleration(101:N+100),10)],0,ones(1,3));
    
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
% figure;
% subplot(1,3,1)
% errorbar(RtoI, Elas, 2*Elas_var)
% hold on
% plot(RtoI, K*ones(length(RtoI),1),'Color','red') 
% subplot(1,3,2)
% errorbar(RtoI, Visc, 2*Visc_var)
% hold on
% plot(RtoI, B*ones(length(RtoI),1),'Color','red') 
% subplot(1,3,3)
% errorbar(RtoI, Iner, 2*Iner_var)
% hold on
% plot(RtoI, I*ones(length(RtoI),1),'Color','red') 

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
%%
[MAG,PHASE,COHE,freq]=get_fresp_IRF(hIntrinsic.dataSet,1/100,4,2);

figure;
subplot(2,1,1)
semilogx(freq,MAG)
subplot(2,1,2)
semilogx(freq,PHASE)

%%
IRF = hIntrinsic.dataSet;
T=1/100;
NLags = 4;
NSides = 2;
[fil_num,fil_den]=besself(5,2*pi*5);  
dis_fil=c2d(tf(fil_num,fil_den),T);

%N=5000;
POSITION=filtfilt(cell2mat(dis_fil.num),cell2mat(dis_fil.den),randn(N,1));
VELOCITY=ddt_real(POSITION)*(1/T);
ACCELERATION=ddt_real(VELOCITY)*(1/T);

IRF_m=irf(irf,'nLags',NLags,'nSides',NSides,'domainIncr',T,...
    'dataSet',IRF,'domainStart',IRF(1));
output=double(nlsim(IRF_m,nldat(POSITION,'domainIncr',T)));


[est,pars]=TV_Bayes(output(100:end-100),[POSITION(100:end-100) VELOCITY(100:end-100) ACCELERATION(100:end-100)]);

INT.K=pars(1);
INT.B=pars(2);
INT.I=pars(3);
INT

[est,pars]=TV_Bayes(tor,[decimate(input_position(101:N+100),10) decimate(velocity(101:N+100),10) decimate(acceleration(101:N+100),10)]);
pars

