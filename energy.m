[sound, fs]=audioread("C:\Users\Amy\Desktop\DSP\Project\noisy_cheese.wav");
t=1/fs*(0:length(sound)-1);
slen=length(sound);
tfinal= t(slen);
frame_length = 0.02; %20ms
frame_shift = 0.01; %10ms

numsamples= slen;
timeps=tfinal/numsamples;
nsamplespf=ceil(0.01/timeps);% of new samples per frame
%entlen= floor(tlen/nsamplespf); %number of frames = length of the entropy array
tsamplespf=floor(0.02/timeps); %num samples pf


%Calculate RMSE
num_frames = floor(numsamples/tsamplespf);
E_max = zeros(1,num_frames);
E_min = zeros(1,num_frames);
VAD = zeros(1,num_frames); %If voice is active/inactive
size(VAD)
parameter = zeros(1,num_frames); %scale parameter for E_min
inactive_count = 0; %counter for inactive frames
entime= (0:tfinal/num_frames:tfinal);
entime= entime(1:length(entime)-1);
threshold= zeros(1,num_frames);
save_energy= zeros(1,num_frames);


for j = 1 :1:num_frames
   %calculating energy of the frame
   current_energy = RMSE(j, num_frames, sound);
   save_energy(1,j)=current_energy;
   
   if(j~=1)
     E_max(1,j)=E_max(1,j-1);
     E_min(1,j)=E_min(1,j-1);
     parameter(1,j) = parameter(1,j-1)*1.0001;
   end
   if (j==1) %if in first frame
     initial = current_energy;
     E_max(1,j)= current_energy;
     E_min(1,j) = current_energy;
     parameter(1,j) = 1.0001; %set initial scale
   elseif (current_energy > E_max(1,j))
     E_max(1,j) = current_energy;
   elseif (current_energy < E_min(1,j))
     if (current_energy <= 0.01)
     E_min(1,j) = initial;
     else
     E_min(1,j) = current_energy;
     end

   end


   %Calculating threshold
   scale = (E_max(1,j) - E_min(1,j))/E_max(1,j);
   threshold(1,j) = ((1-scale)*E_max(1,j)) + (scale*E_min(1,j));
   if(current_energy > threshold(1,j))
     VAD(1,j) = 1;
     inactive_count = 0;
   elseif (inactive_count == 4)
     VAD(1,j) = 0;
   else
     VAD(1,j) = 1;
     inactive_count = inactive_count + 1;
   end
   E_min(1,j) = E_min(1,j)*parameter(1,j); %%not letting me do E_min(j) = E_min(j1)*parameter(j);

end


VAD;
figure();

% plot(t,sound);
plot(entime, save_energy, 'b');
% size(entime)
% size(E_max)
hold on;
plot(entime,E_max);
hold on;
plot(entime,E_min);
hold on;
plot(entime, VAD);
hold on;
plot(entime, threshold, 'color','g')
%plot(sound);
%ylabel('low frequency');
%hold on;
% figure


function energy = RMSE(j,N, signal)
   start = ((j-1)*N)+1;
   final = j*N;
   energy = 0;
   for i = start:1:final
     energy = energy + (signal(i).^2);
   end
   energy = energy/N;
   energy = sqrt(energy);
end