function [seq] = gen_prbs(SI,amp,length,fs)
% randn('state',0);
% rand('state',0); 

% generate a pseudo-random binary sequence
dLength_s = length;
dAmplitude = amp;
dInterval_ms = SI;
nDARate_Hz = fs;

dIncrement = 1000 / nDARate_Hz;
    
nSwitch = round(dLength_s * 1000 / dInterval_ms)';
% nSwitch = round(nSwitch);

Stimulus = ( rand ( nSwitch, 1 ) >= .5 );
Stimulus = Stimulus * dAmplitude / range( Stimulus ); 

% Stimulus = ( Stimulus - mean( Stimulus ) ); %include for no random delay

np = dInterval_ms / dIncrement;

SequenceStimulus = [];

for i=1:nSwitch, %length( Stimulus ),
    SequenceStimulus = [ SequenceStimulus zeros( 1, np ) + Stimulus( i ) ];
end

%% HIG: Added random delay
dDelay = round(rand*np);
dIC = (rand(1,1) >= .5)*dAmplitude;
SequenceStimulus = [ones(1,dDelay)*dIC SequenceStimulus(1,1:end-dDelay)]; 
SequenceStimulus = SequenceStimulus - mean(SequenceStimulus); 
%%

seq = SequenceStimulus';