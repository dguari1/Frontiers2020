T = 1/1000;
delay = 40/1000;

N = 30000;

INPUT_POS=load('C:\Users\GUARIND\Dropbox (Personal)\PhD\Paper_JointStiff_Movement\Code\position_input.mat');
input_position=INPUT_POS.position_input;



PRBS = gen_prbs(250,0.03*2,N*T,1/T);
fc=3;% cut off frequency
fn=1000/2; %nyquivst frequency = sample frequency/2;
order = 2; %6th order filter, high pass
[b_lp, a_lp]=butter(order,(fc/fn),'low');
input_position = filtfilt(b_lp,a_lp, PRBS);

velocity=ddt_real(input_position)*(1/T);
acceleration=ddt_real(velocity)*(1/T);

figure;subplot(2,1,1);plot(input_position);subplot(2,1,2);plot(velocity)
%%
NOISE=load('C:\Users\GUARIND\Dropbox (Personal)\PhD\Paper_JointStiff_Movement\Code\Noise.mat');
out_noise=NOISE.Noise;


K = 40;
B = 1.5;
I = 0.02;

Intrinsic_torque = K*input_position(1:N) + B*velocity(1:N) + I*acceleration(1:N);
Intrinsic_torque = Intrinsic_torque(1:N);


tau=ceil(delay/T); %discrete delay
del_velocity=[zeros(tau,1);velocity(1:(N)-tau)];  %delayed velocity
nl_del_velocity = del_velocity;
nl_del_velocity(nl_del_velocity<0.25)=0;

gainn= 0:0.5:20;

G = -14.5;
w = 16.1;
z = 0.65;

%to split data 
[pks,locs]=findpeaks(-(velocity(1:N)),'MinPeakHeight',1,'MinPeakDistance',100);

len = 40:10:240;

for k =1 : length(len)%length(gainn)
    G = -3;%-gainn(k);
    Reflex_torque = lsim(tf(G*w*w,[1,2*z*w,w*w]), nl_del_velocity, (0:N-1)*T);

    
    Reflex_torque = Reflex_torque(1:N);


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

    %split data in segments of +40 /\ -40 seg
    clear ens_pos ens_vel ens_acc ens_tor ens_int
    S=len(k);%40;
    for i = 1:length(locs)

        ens_pos(:,i) = decimate(input_position(locs(i)-S:locs(i)+S),10);

        ens_vel(:,i) = decimate(velocity(locs(i)-S:locs(i)+S),10);

        ens_acc(:,i) = decimate(acceleration(locs(i)-S:locs(i)+S),10);

        ens_tor(:,i) = decimate(Total_torque(locs(i)-S:locs(i)+S),10);
        
        ens_int(:,i) = decimate(Intrinsic_torque(locs(i)-S:locs(i)+S),10);

    end
    
    [Params_int,Ps_int,est] = my_rivbjmiso(vec(ens_tor),[vec(ens_pos) vec(ens_vel) vec(ens_acc)],0,ones(1,3));

    Elas_riv(k,1) = Params_int(1);
    Visc_riv(k,1) = Params_int(2);
    Iner_riv(k,1) = Params_int(3);
    
    Elas_var_riv(k,1) = sqrt(Ps_int(1,1));
    Visc_var_riv(k,1) = sqrt(Ps_int(2,2));
    Iner_var_riv(k,1) = sqrt(Ps_int(3,3));
    
    m = length(ens_pos(:,1));
    ens_est  = vec2mat(est,m)';
    clear VAFy_tor VAFy_int
    for p = 1 :length(locs)
        VAFy_tor(p,1) = VAFnl(ens_tor(:,p),ens_est(:,p));
        VAFy_int(p,1) = VAFnl(ens_int(:,p),ens_est(:,p));
    end
    
     VAF_total_mean(k,1) = mean(VAFy_tor);
     VAF_total_std(k,1) = std(VAFy_tor);
     VAF_int_mean(k,1) = mean(VAFy_int);
     VAF_int_std(k,1) = std(VAFy_int);  

end

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
errorbar(len, Elas_riv, 2*Elas_var_riv)
hold on
plot(len, K*ones(length(len),1),'Color','red') 
subplot(1,3,2)
errorbar(len, Visc_riv, 2*Visc_var_riv)
hold on
plot(len, B*ones(length(len),1),'Color','red') 
subplot(1,3,3)
errorbar(len, Iner_riv, 2*Iner_var_riv)
hold on
plot(len, I*ones(length(len),1),'Color','red') 
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

