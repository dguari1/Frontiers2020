
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

clear RtoI Elast_est Visc_est Iner_est VAF_total VAF_int IRF

max_iteration = 10;
id_tolerance = 0.05;
mReflex=nlbl;
set(mReflex,'idMethod','sls', 'displayFlag',false,'nIterMax',20,'nLagLE',80, 'maxOrderNLE',6);
hIntrinsic = irf;
hIntrinsic= set(hIntrinsic,'nLags',4,'nSides',2);


for m = 1 : 100 %100 simulation trials 
    for k =1 : length(gainn)
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
        G = -gainn(k);
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


        %estimation
        
        tqI_temp = zeros(N/10,1);
        tqR_temp = zeros(N/10,1);

        VAFbest = -1000;
        
        pos = decimate(input_position,10);
        vel = decimate(velocity,10);
        tor = decimate(Torque,10);
        for iter = 1 : max_iteration
            tqI_temp = tor - tqR_temp;

%             zp=cat(2,nldat(pos,'domainIncr',1/100),nldat(tqI_temp,'domainIncr',1/100));
%             hIntrinsic = nlident(hIntrinsic,zp);
% 
%             tqI_temp= double(nlsim(hIntrinsic,nldat(pos,'domainIncr',1/100)));
%             
            Params_int = [decimate(input_position(1:N),10) decimate(velocity(1:N),10) decimate(acceleration(1:N),10)]\tqI_temp;
            tqI_temp = [decimate(input_position(1:N),10) decimate(velocity(1:N),10) decimate(acceleration(1:N),10)]*Params_int;
            
            if k >1
                tqR_temp = tor-tqI_temp;
                zv=cat(2,nldat(vel,'domainIncr',1/100),nldat(tqR_temp,'domainIncr',1/100));
                mReflex=nlident(mReflex, zv);
                tqR_temp= double(nlsim(mReflex,nldat(vel,'domainIncr',1/100)));

                tqT_temp = tqI_temp + tqR_temp;

                VAFT = VAFnl(tor,tqT_temp);

                if (VAFT-VAFbest) < id_tolerance
                    break
                else
                    VAFbest=VAFT;
                end
            else
                tqT_temp = tqI_temp;
                break
            end
        end
% 
%         zp = cat(2,nldat(input_position,'domainIncr',1/1000), nldat(Torque,'domainIncr',1/1000));
%         [hIntrinsic, mReflex, nlmStiff, VAFT, VAFi, VAFr, tqI, tqR, tqT]= PC_stiffnessID (zp, 'reflex_irf_len',300,'reflex_id_method', 'sls');
%         
%        est =  double(tqI_temp);
        
%        Params_int = [decimate(input_position(1:N),10) decimate(velocity(1:N),10) decimate(acceleration(1:N),10)]\est;
        %est = [decimate(input_position(1:N),10) decimate(velocity(1:N),10) decimate(acceleration(1:N),10)]*Params_int;
        Elast_est(k,m) = Params_int(1);
        Visc_est(k,m) = Params_int(2);
        Iner_est(k,m) = Params_int(3);

%         Elas_var_riv(k,1) = sqrt(var(Total_noise))*sqrt(Ps_int(1,1));
%         Visc_var_riv(k,1) = sqrt(var(Total_noise))*sqrt(Ps_int(2,2));
%         Iner_var_riv(k,1) = sqrt(var(Total_noise))*sqrt(Ps_int(3,3));
% 

        VAF_total(k,m) = VAFnl(decimate(Total_torque,10),double(tqT_temp));
        VAF_int(k,m) = VAFnl(decimate(Intrinsic_torque,10),tqI_temp);


    end
    
end
%%
Mean_Elast = mean(Elast_est');
Std_Elast = std(Elast_est');

Mean_Visc = mean(Visc_est');
Std_Visc  = std(Visc_est');

Mean_Iner = mean(Iner_est');
Std_Iner = std(Iner_est');

figure;
subplot(1,3,1)
plot(mean(RtoI'), K*ones(length(gainn),1),'Color','red') 
hold on
errorbar(mean(RtoI'), Mean_Elast', 1.5*Std_Elast')
title('Elasticity')
xlabel('ratio_{RtoI}')
ylabel('K (Nm/rad)')

subplot(1,3,2)
plot(mean(RtoI'), B*ones(length(gainn),1),'Color','red') 
hold on
errorbar(mean(RtoI'), Mean_Visc, 1.5*Std_Visc)
title('Viscosity')
xlabel('ratio_{RtoI}')
ylabel('B (Nm/rad/s)')


subplot(1,3,3)
plot(mean(RtoI'), I*ones(length(gainn),1),'Color','red') 
hold on
errorbar(mean(RtoI'), Mean_Iner', 1.5*Std_Iner')
title('Inertia')
xlabel('ratio_{RtoI}')
ylabel('I (Nm/rad/s^2)')


Mean_VAFtotal = mean(VAF_total');
Std_VAFtotal = std(VAF_total');

Mean_VAFint = mean(VAF_int');
Std_VAFint  = std(VAF_int');

figure
subplot(1,2,1)
errorbar(mean(RtoI'), Mean_VAFtotal', 2*Std_VAFtotal')
title('%VAF  - Overall Torque')
xlabel('ratio_{RtoI}')
ylabel('%VAF')
subplot(1,2,2)
errorbar(mean(RtoI'), Mean_VAFint', 2*Std_VAFint')
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

