[sound, fs]=audioread("C:\Users\Amy\Desktop\DSP\Project\dirty\this_moment_white.wav");
[origsound, fs1]=audioread("C:\Users\Amy\Desktop\DSP\Project\clean\this_moment.wav");
sound= sound(:,1);
origsound= origsound(:,1);
sound = awgn(sound,2);
t=1/fs*(0:length(sound)-1);
to=1/fs1*(0:length(origsound)-1);
slen=length(sound);
tfinal= t(slen);
frame_length = 0.02; %20ms
frame_shift = 0.01; %10ms

dnumsamples= slen;
timeps=tfinal/numsamples;
nsamplespf=floor(0.01/timeps);% of new samples per frame
tsamplespf=floor(0.02/timeps); %num samples pf

W = 20*dbaux(8); %mother wavelet DB8 length =16
hscales = [5 6 7 8 9 10 15];
lscales=[40 60 90 110 120 150 300 500];

LCOEFF = WaveConv(sound, lscales, W,t);

[LCE, Ltime] = GetEntropyVector(LCOEFF,tsamplespf,nsamplespf, timeps);
rows= size(LCOEFF);
netl= zeros(1, length(LCOEFF));
for i = 1 :1:rows(1)
  netl= netl+ LCOEFF(i,:);
end
ltimeps= tfinal/length(LCOEFF);


[realent, rt]=GetEntropyVector(netl,tsamplespf,nsamplespf, ltimeps);


realent= realent.^10;
realent= smooth(realent,0.01,'moving');
realent=realent.^4;


%_____________________________________________________________________
%CALC THRESHOLD
%%%find frames
slen=length(realent);
tfinal= rt(slen);
frame_length = 0.02; %20ms
frame_shift = 0.01; %10ms


numsamples= slen;
timeps=tfinal/numsamples;
nsamplespf=ceil(0.01/timeps);% of new samples per frame
%entlen= floor(tlen/nsamplespf); %number of frames = length of the entropy array
tsamplespf=floor(0.02/timeps); %num samples pf

num_frames = numsamples;
E_max = zeros(1,num_frames);
E_min = zeros(1,num_frames);
maxi=zeros(1,num_frames);
mini=zeros(1,num_frames);
VAD = zeros(1,num_frames); %If voice is active/inactive
parameter = zeros(1,num_frames); %scale parameter for E_min
inactive_count = 0; %counter for inactive frames
threshold= zeros(1,num_frames);
siz=size(realent);

for j = 1 :1:length(realent)

   current_ent = realent(j);
   %current_ent=transpose(current_ent);
   if(j~=1)
     E_max(1,j)=E_max(1,j-1);
     E_min(1,j)=E_min(1,j-1);
     parameter(1,j) = parameter(1,j-1);
   end
   if (j<=20) %if in first frame
     initial = current_ent;
     E_max(1,j)= 0.1;
     E_min(1,j) = 0.8;
     parameter(1,j) = 1.0001; %set initial scale
   elseif (current_ent > E_max(1,j))
     E_max(1,j) = current_ent;
   elseif (current_ent < E_min(1,j))
     if (current_ent == 0)
     E_min(1,j) = initial;
     else
     E_min(1,j) = current_ent;
     end

   end

   %Calculating threshold
   scale = (E_max(1,j) - E_min(1,j))/E_max(1,j);
   threshold(1,j) = ((1-scale)*E_max(1,j)) + (scale*E_min(1,j));
   if(current_ent < threshold(1,j))
     VAD(1,j) = 1;
     inactive_count = 0;
   elseif (inactive_count == 6*tsamplespf)
     VAD(1,j) = 0;
   else
     VAD(1,j) = 1;
     inactive_count = inactive_count + 1;
   end
   E_min(1,j) = E_min(1,j)*parameter(1,j);
   if j<5
     VAD(1,j)=0;
   end
end

% maxi
% mini
entime= (0:tfinal/num_frames:tfinal);
entime= entime(1:length(entime)-1);
figure();
plot(to, origsound);
%ylabel('low frequency');
hold on;
plot(rt,realent,'r');
hold on;
plot(rt, VAD, 'g');
% hold on;
% plot(rt, threshold, 'm');
% % % plot(entime, E_min);
% hold on;
% plot(rt, E_min, 'c');


function CD = WaveConv(sound, sc, mw,t)
   CD = zeros(length(sc), length(sound));
   cds=size(CD);
   for m = 1:1:length(sc)
       s= sc(m);
       %dbaux is real valued so conj(dbaux)= dbaux
       x = 1:1:16;
       v=mw;
      %Define the query points to be a finer sampling over the range of x.
       xq = 1:1/s:16;
      % Interpolate the function at the query points and plot the result.
       tmw= interp1(x,v,xq); % time scaled : mv(t/s)
       f = 1/sqrt(s) * tmw ; %reference wavelet = shifted, scaled version of mother
      wavelet s^(-1/2)* mv(t/s) page 12 doc
       CF = conv(sound,f);
       CF = CF(1:length(sound)+1);
       lcf = length(CF);
       %CD(m,:) = diff(CF)./diff(t);
   end
end


function [e, te]= GetEntropyVector(i,fl, fs, timeps) %sound i, fl=totalsamplespf,fs=newsamplespf, timepersample
   len = length(i);
   sp = 1; %start point
   ep = sp + fl ; %end point = start point + frame length
   count=1;
   %initialise e and te to be too big
   e=zeros(1,len);
   te=zeros(1,len);
   while(ep <= len)
     t=0:1:fl;
     t=t+sp-1;
     t=t.*timeps;
     inp=i(sp:ep); %find frame
     [ent,time] = pentropy(inp,t);%calculate entropy of frame
     ent=transpose(ent);
     totsize=length(ent);
     %%% now we combine the entropy frames
     %find the overlap length on each side of entropy frames
     % overlap total =10ms. overlap on each side is 5ms
     %overlap is 5ms/20ms= 1/4 of a frame
     overlap=floor(totsize/4);
     if count~=1
       %read last 5ms of entropy from the entropy vector
       last5= e(count-overlap-1: count-1);
       %read first 5ms from entropy calculation
       first5=ent(1:1+overlap);
       %find the average of the entropy calculated and the entropy
       %saved which overlap each other
       avg= (first5+last5)./2;
       %save average to the entropy vector
       e(count-overlap-1:count-1)=avg;
       %save remaining calculated entropy to the entropy vector
       e(count:count+totsize-overlap-1)= ent(overlap+1:totsize);
       %save time to time vector
       te(count:count+totsize-overlap-1)= time(overlap+1:totsize);
     else % if the count ==1, save full entropy and time as initial condition
       e(1: totsize)=ent;
       te(1:totsize)=time;
     end
     % increase sp and ep
     sp = sp + fs;
     ep = sp + fl ;
     %move the count value to the end of the frame
     count=count+totsize;
   end
   % fix the sizes of e and te by cropping off the spare zeros
   e= e(1:count-1);
   te=te(1:count-1);

end