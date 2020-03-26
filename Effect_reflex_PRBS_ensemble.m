
%simulation parameters 
T = 1/1000;
N = 60000;

%model parameters 
%intrinsic
K = 30;
B = 0.8;
I = 0.012;

%reflex 
gainn= 0:0.5:10;
w = 15.6;
z = 0.62;
ref_threshold = 0.5;
delay = 40/1000; 

clear RtoI Elast_est Visc_est Iner_est VAF_total VAF_int
        
len = 40:10:240;
for m = 1 : 100 %100 simulation trials 
    for k =1 : length(len)
        clear PRBS input_position velocity acceleration Intrinsic_torque Reflex_torque Total_torque Torque est
        %compute inputs 
        PRBS = gen_prbs(200,0.03*2,N*T,1/T);
        fc=30;% cut off frequency
        fn=1000/2; %nyquivst frequency = sample frequency/2;
        order = 2; %6th order filter, high pass
        [b_lp, a_lp]=butter(order,(fc/fn),'low');
        input_position = filtfilt(b_lp,a_lp, PRBS);
        %figure;subplot(2,1,1);plot((0:T:(N-1)*T),input_position);subplot(2,1,2);plot((0:T:(N-1)*T),gradient(input_position,T))
        %figure;plot(spect(nldat(PRBS-mean(PRBS),'domainIncr',T)))
        velocity=gradient(input_position,T);
        acceleration=gradient(velocity,T);

        %intrnisic torque 
        Intrinsic_torque = K*input_position(1:N) + B*velocity(1:N) + I*acceleration(1:N);
        Intrinsic_torque = Intrinsic_torque(1:N);

        %compute delayed velocity 
        tau=ceil(delay/T); %discrete delay
        del_velocity=[zeros(tau,1);velocity(1:(N)-tau)];  %delayed velocity
        nl_del_velocity = del_velocity;
        nl_del_velocity(nl_del_velocity<ref_threshold)=0;

    
        %reflex torque 
        G = -5.50;
        Reflex_torque = lsim(tf(G*w*w,[1,2*z*w,w*w]), nl_del_velocity, (0:N-1)*T);    
        Reflex_torque = Reflex_torque(1:N);

        %computer ratio RtoI
        RtoI(k,m) = var(Reflex_torque)/var(Intrinsic_torque);

        Total_torque = Intrinsic_torque+Reflex_torque;


        %noise
        Noise = randn(N,1);
        power_noise = sum(Noise.^2);
        power_signal = sum((Total_torque).^2);
        snr=15;
        Total_noise = Noise*sqrt((power_signal/(10^(snr/10)))/power_noise);

        Torque = Total_torque + Total_noise;

        [pks,locs]=findpeaks(-(velocity(1:N)),'MinPeakHeight',1,'MinPeakDistance',100);
        clear ens_pos ens_vel ens_acc ens_tor ens_int
        S=len(k);%40;
        for i = 2:length(locs)-1
            ens_pos(:,i) = decimate(input_position(locs(i)-S:locs(i)+S),10);

            ens_vel(:,i) = decimate(velocity(locs(i)-S:locs(i)+S),10);

            ens_acc(:,i) = decimate(acceleration(locs(i)-S:locs(i)+S),10);

            ens_tor(:,i) = decimate(Total_torque(locs(i)-S:locs(i)+S),10);
        end
        
        Params_int = [vec(ens_pos) vec(ens_vel) vec(ens_acc)]\vec(ens_tor);
        
        %Params_int = [decimate(input_position(1:N),10) decimate(velocity(1:N),10) decimate(acceleration(1:N),10)]\decimate(Torque,10);
        est = [decimate(input_position(1:N),10) decimate(velocity(1:N),10) decimate(acceleration(1:N),10)]*Params_int;
        Elast_est(k,m) = Params_int(1);
        Visc_est(k,m) = Params_int(2);
        Iner_est(k,m) = Params_int(3);

%         Elas_var_riv(k,1) = sqrt(var(Total_noise))*sqrt(Ps_int(1,1));
%         Visc_var_riv(k,1) = sqrt(var(Total_noise))*sqrt(Ps_int(2,2));
%         Iner_var_riv(k,1) = sqrt(var(Total_noise))*sqrt(Ps_int(3,3));
% 

        VAF_total(k,m) = VAFnl(decimate(Total_torque,10),est);
        VAF_int(k,m) = VAFnl(decimate(Intrinsic_torque,10),est);


    end
    
end

Mean_Elast = mean(Elast_est');
Std_Elast = std(Elast_est');

Mean_Visc = mean(Visc_est');
Std_Visc  = std(Visc_est');

Mean_Iner = mean(Iner_est');
Std_Iner = std(Iner_est');

figure;
subplot(1,3,1)
plot(len, K*ones(length(len),1),'Color','red') 
hold on
errorbar(len, Mean_Elast', 2*Std_Elast')
title('Elasticity')
xlabel('ratio_{RtoI}')
ylabel('K (Nm/rad)')

subplot(1,3,2)
plot(len, B*ones(length(len),1),'Color','red') 
hold on
errorbar(len, Mean_Visc, 2*Std_Visc)
title('Viscosity')
xlabel('ratio_{RtoI}')
ylabel('B (Nm/rad/s)')


subplot(1,3,3)
plot(len, I*ones(length(len),1),'Color','red') 
hold on
errorbar(len, Mean_Iner', 2*Std_Iner')
title('Inertia')
xlabel('ratio_{RtoI}')
ylabel('I (Nm/rad/s^2)')


Mean_VAFtotal = mean(VAF_total');
Std_VAFtotal = std(VAF_total');

Mean_VAFint = mean(VAF_int');
Std_VAFint  = std(VAF_int');

figure
subplot(1,2,1)
errorbar(len, Mean_VAFtotal', 2*Std_VAFtotal')
title('%VAF  - Overall Torque')
xlabel('ratio_{RtoI}')
ylabel('%VAF')
subplot(1,2,2)
errorbar(len, Mean_VAFint', 2*Std_VAFint')
title('%VAF  - Intrinsic Torque')
xlabel('ratio_{RtoI}')
ylabel('%VAF')
% %%
% [MAG,PHASE,COHE,freq]=get_fresp_IRF(hIntrinsic.dataSet,1/100,4,2);
% 
% figure;
% subplot(2,1,1)
% semilogx(freq,MAG)
% subplot(2,1,2)
% semilogx(freq,PHASE)
% 
% %%
% IRF = hIntrinsic.dataSet;
% T=1/100;
% NLags = 4;
% NSides = 2;
% [fil_num,fil_den]=besself(5,2*pi*5);  
% dis_fil=c2d(tf(fil_num,fil_den),T);
% 
% %N=5000;
% POSITION=filtfilt(cell2mat(dis_fil.num),cell2mat(dis_fil.den),randn(N,1));
% VELOCITY=ddt_real(POSITION)*(1/T);
% ACCELERATION=ddt_real(VELOCITY)*(1/T);
% 
% IRF_m=irf(irf,'nLags',NLags,'nSides',NSides,'domainIncr',T,...
%     'dataSet',IRF,'domainStart',IRF(1));
% output=double(nlsim(IRF_m,nldat(POSITION,'domainIncr',T)));
% 
% 
% [est,pars]=TV_Bayes(output(100:end-100),[POSITION(100:end-100) VELOCITY(100:end-100) ACCELERATION(100:end-100)]);
% 
% INT.K=pars(1);
% INT.B=pars(2);
% INT.I=pars(3);
% INT
% 
% [est,pars]=TV_Bayes(tor,[decimate(input_position(101:N+100),10) decimate(velocity(101:N+100),10) decimate(acceleration(101:N+100),10)]);
% pars
% 
