%read data
data  = flb2mat('C:\Users\guarind\Dropbox (NRP)\Intrinsic-reflex EMG\MR_040412.flb','read_case',17);

%sampling rate
T=1/1000;
%take the what data will be used
POSITION = data.Data(:,1) - mean(data.Data(:,1));
VELOCITY = ddt_real(POSITION,T);
ACCELERATION = ddt_real(VELOCITY,T);
TORQUE = data.Data(:,2) - mean(data.Data(:,2));
%decimate data to remove high-freq noise
dr= 10;
position = decimate(POSITION,dr);
velocity = decimate(VELOCITY,dr);
acceleration = decimate(ACCELERATION,dr);
torque = decimate(TORQUE,dr);

%separate into identification (id) and validation (val) data sets
id_pos = position(1:5000);
id_vel = velocity(1:5000);
id_acc = acceleration(1:5000);
id_tor = torque(1:5000);

val_pos = position(5001:end);
val_vel = velocity(5001:end);
val_acc = acceleration(5001:end);
val_tor = torque(5001:end);

%identify intrinsic model (IBK) parameters
[pred_id,pars,noise_var,lk,COV]=TV_Bayes(id_tor,[id_pos, id_vel, id_acc]);

INT.K=pars(1);
INT.B=pars(2);
INT.I=pars(3);

%compute identification VAF
VAF_id = VAFnl(id_tor, pred_id);

%validate the model
pred_val = INT.K*val_pos + INT.B*val_vel + INT.I*val_acc;
VAF_val = VAFnl(val_tor, pred_val);

[-INT.K, -INT.B, -INT.I]