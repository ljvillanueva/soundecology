[sig,fs] = wavread('ROMO_20100518_020000.wav');
[bL1,f,tL1] = specgram(sig(:,1),1024,fs); % left side first (L)
dbL1=20*log10(abs(bL1));
[bR1,f,tR1] = specgram(sig(:,2),1024,fs);  % right side second (R)
dbR1=20*log10(abs(bR1));
clear sig
dbL1_avg=mean(dbL1');                 % avg over time for L
dbR1_avg=mean(dbR1');                  % avg over time for R
db1_LR_avg=[dbL1_avg;dbR1_avg];

%next, average the L and R together
db1_avg=mean(db1_LR_avg,1);
save ROMO_20100518_020000.dat db1_avg -ascii
load ROMO_20100518_020000.dat
clear db*
clear b*
clear t*

%Hz=[180:38.64:20000];
%plot(Hz,ROMO_20100518_020000,'b')

% AREA calcultion for full spectrum
%freq= diff(Hz);
%y1=ROMO_20100518_020000 - min(ROMO_20100518_020000);           
%y1= y1(1:length(y1)-1);
%ROMO_20100518_020000_area= sum(y1.*freq);             

%clear y*
%save ROMO_20100518_020000_area.dat ROMO_20100518_020000_area -ascii


%CUT THE SPECTRUM DOWN TO 180Hz through 1500Hz and PLOT it

ROMO_20100518_020000_short=ROMO_20100518_020000(:,35:217);  % selects only values from 180 to 8500 Hz
save ROMO_20100518_020000_short.dat ROMO_20100518_020000_short -ascii
Hz_short=[1500:38.35:8500];    %makes the Hz files to plot the above against
figure
plot(Hz_short,ROMO_20100518_020000_short)

%AREA INDEX for SHORT spectrum
freq_short= diff(Hz_short);
y1=ROMO_20100518_020000_short - min(ROMO_20100518_020000_short);           
y1= y1(1:length(y1)-1);
ROMO_20100518_020000_short_area= sum(y1.*freq_short);  


%SUM of all dB at all frequencies between 1996 Hz and 
