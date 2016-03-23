clc;
close all;
clear all;
%% If noise used at the channel then BER curve is good but offset estimated is not good.
q=1;w=1;e=1;r=1;
%% Characterisng channel
fs=100;
ts=1/fs; 
%% Data
M=4;
for N=[2^8  2^9  2^10  2^11]
Nfft=N;
Nps=4; % Pilot spacing
No_pilot=N/Nps;
%% Modulation
data=randi([0 M-1],1,N-No_pilot);
mod_data=qammod(data,M);ep_new=zeros(1,30);ep=ones(1,30);ep=0.15*(ep);error_offset=zeros(1,30);

%% Pilot Insertion
Nps=4; % Pilot spacing
No_pilot=N/Nps;
Xp = 2*(randn(1,No_pilot)>0)-1;% Pilot seq generation
ip = 0; pilot_loc = [];
for k=1:N
if mod(k,Nps)==1
X(k)=Xp(floor(k/Nps)+1); 
pilot_loc=[pilot_loc k]; 
ip = ip+1;
else X(k) = mod_data(k-ip);
end
end

Nfft=N; 
x = ifft(X,Nfft); % IFFT of pilot inserted data
%% Add Cyclic Perfix
Ng=N/4;
xt = [x(Nfft-Ng+1:Nfft) x];
 x=xt;
 %% Providing frequency offset
 ep1=0.15;
 mm=(0:(length(x))-1);
 offset = exp(1i*2*pi*ep1/(length(x)).*mm);
 x=x.*offset;
 
 %% Passing sig through channel

 h = [(randn+j*randn) (randn+j*randn)/2 (randn+j*randn)/3 (randn+j*randn)/4]; % A (4-tap) channel
H = fft(h,Nfft); 
ch_length=length(h); % True channel and its length
H_power_dB = 10*log10(abs(H.*conj(H))); % True channel power in dB
y_channel = conv(x,h); % Channel path (convolution)
%% Noise
    for snr_db=1:30
    sig_power=bandpower(data);
     snr_l=10.^(snr_db/10);
     var=sig_power/snr_l;
     noise=(1/sqrt(2)*(randn(1,length(y_channel))+j*randn(1,length(y_channel))))*var;
     noisy_sig=y_channel+noise;
     noisy_sig=awgn(y_channel,snr_db,'measured');
%% Remove Cyclic Prefix
 Nofdm=N+Ng;
 Y = noisy_sig(Ng+1:Nofdm);
 Y=fft(Y,Nfft);

%% Estimation of Channel
m=1;
H_est = LS_CE(Y,Xp,pilot_loc,Nfft,Nps,m);
H_est_power_dB = 10*log10(abs(H_est.*conj(H_est)));
h_est = ifft(H_est);
h_DFT = h_est(1:ch_length);
H_DFT = fft(h_DFT,Nfft); % DFT-based channel estimation
H_DFT_power_dB = 10*log10(abs(H_DFT.*conj(H_DFT)));
Y_eq=Y./H_est; 
%% Removing Pilot
ip = 0;
for k=1:Nfft
if mod(k,Nps)==1, 
    ip=ip+1;
else
    Data_extracted(k-ip)=Y_eq(k);
end
end
 %% Calculate offset
      mm=(0:(N-No_pilot)-1);
      ep_new(snr_db)=new_extended_kalmaan((Data_extracted),(data),var);
       offset_cal = exp(-1i*2*pi*ep_new(snr_db)/(N-No_pilot).*mm);
      de_offset_data=Data_extracted.*offset_cal;
        
     
%% Demodulation
 dem_data=qamdemod(de_offset_data,M);
 
 %% Error cal
 BER(snr_db)=symerr(data,dem_data)/N;
 error_offset(snr_db)=((ep(snr_db)-ep_new(snr_db))/ep(snr_db))*100;

 
 %% Plotting
 if snr_db==30
     figure(1);
    %% ploting channel response
     subplot(2,2,q)
     method='LS-linear';
     plot(H_power_dB,'r'); 
     hold on
    plot(H_DFT_power_dB,'b--');
    legend('True Channel',[method  'with DFT']);
    ylabel('|H(f)|')
    str=sprintf('Plot with SNR = %d and no of samples = %d', 30, N);
    title(str)
    hold off
    %% ploting data vs rcvd_data
    figure(2);
    subplot(2,2,w)
    stem(data(1:20),'r:+')
    hold on
    stem(dem_data(1:20),'b')
    legend('transmitted data','received data')
    str=sprintf(' Tx and Rx and no of samples = %d', N);
    title(str)
    hold off
    w=w+1;
    q=q+1;
 end
    end
    n=1:30;
%     snr_l=10.^(n./10);
    BER_th=(1/2)*erfc(sqrt(snr_l)); 
    figure(3);
    %% plotting BER curve
    subplot(2,2,e)
     semilogy(n,BER,'*r')
    hold on
%     semilogy(n,BER_th,'k*');
    xlabel('SNR dB')
    ylabel('BER')
    str=sprintf(' BER curve with no of samples = %d', N);
    title(str)
     legend('Simulation');
    e=e+1;
     figure(4);
     subplot(2,2,r)
     semilogy(n,error_offset,'k*')
     hold on
     semilogy(n,BER_th,'k*');
     xlabel('SNR dB')
     ylabel('% error in Offset')
     str=sprintf(' per error vs snr-dB and no of samples = %d and offset = %d', N,max(ep));
     title(str)
     
      legend('% error in offset');
     r=r+1;
     
    
end
   

