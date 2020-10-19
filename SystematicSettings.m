function [SNR,order,TC,PDT,CNR]=SystematicSettings(CNR,PDT,pixelsize,objectsize,Pcarrier,Pprobe,Csignal,Cshot,Cthermal,choice,order,ordermax,w01x,w01y,w02x,w02y,zR1x,zR2x,zR1y,zR2y)
% Written by Bart Fokker, Vrije Universiteit Amsterdam, 10 January 2019
%edited by: Liron Zada. Last edit: 14 may 2020.
% Purpose: calculate the optimal settings for our SRS setup.

% Initialize parameters
orders=(order:ordermax)';
OptimalTC=zeros(length(orders),1);
CNRproportional=zeros(length(orders),1);
response=zeros(length(orders),1);
bestoutput=cell(length(orders),1);

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% L.Z.:
% First, define TC range for each filter order
% Secondly,find the optimal "pre-TC" with LIAOptimalSettings
% Thirdly, find the real TC and optimal CNR or TC and PDT. 
% Fourthly, Plot signal response 
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% Run the simulation for each filter order
for i=1:length(orders)
    
    % Set the available TCs for each filter order
    % L.Z: These TC ranges were selected to be proportional to 100.85 us PDT
    switch orders(i)
        case 1
            TC=(20*10^-6:10^-6:160*10^-6)';
        case 2
            TC=(10*10^-6:10^-6:100*10^-6)';
        case 3
            TC=(5*10^-6:10^-6:100*10^-6)';
        case 4
            TC=(5*10^-6:10^-6:50*10^-6)';
        case 5
            TC=(5*10^-6:10^-6:50*10^-6)';
        case 6
            TC=(5*10^-6:10^-6:40*10^-6)';
        case 7
            TC=(5*10^-6:10^-6:40*10^-6)';
        otherwise
            TC=(1*10^-6:10^-6:30*10^-6)';
    end

    % Calculate optimal pre-TC settings
    if i~=length(orders)
        [OptimalTC(i),CNRproportional(i),response(i),bestoutput{i}]=LIAOptimalSettings(objectsize,choice,orders(i),TC,w01x,w01y,w02x,w02y,zR1x,zR2x,zR1y,zR2y);
    else
        [OptimalTC(i),CNRproportional(i),response(i),bestoutput{i},sample,InputSample,samplesperpixel,numberofpixels,beamprecision]=LIAOptimalSettings(objectsize,choice,orders(i),TC,w01x,w01y,w02x,w02y,zR1x,zR2x,zR1y,zR2y);
    end
end

% Extract optimal TC and filter order
[~,Index] = max(CNRproportional);
BestTC = OptimalTC(Index);
order = orders(Index);
ChosenResponse = response(Index); % For corrected output values

% Convert pre-TC and input parameters to TC, PDT and CNR
if isempty(CNR) %L.Z. CNR optimization. 
    TC=BestTC*PDT*objectsize/(100.85*10^-6*pixelsize);
    SNR=Csignal*Pprobe*Pcarrier*sqrt(4*sqrt(pi)*TC*gamma(order)/(gamma(order-1/2)*(Cthermal^2+Pprobe*Cshot^2)));
    
    % Corrected output values; L.Z.ChosenResponse converts SNR to CNR
    CNR=SNR*ChosenResponse/sqrt(2); %L.Z. sqrt(2) is added to account for the double noise contribution. 
    
else %L.Z. Imaging time optimization. 
    
     % When Contrast is specified (here being the initial CNR parameter)
     SNR=CNR*sqrt(2)/ChosenResponse;
     TC=(SNR/(Csignal*Pprobe*Pcarrier))^2*(gamma(order-1/2)*(Cthermal^2+Pprobe*Cshot^2))/(4*sqrt(pi)*gamma(order))/ChosenResponse^2;
     PDT=TC*100.85*10^-6*pixelsize/(BestTC*objectsize);
    
end


% Plot sample, input, and optimal output
figure

% Plot object
sampleX=(1/round(objectsize/beamprecision):1/round(objectsize/beamprecision):numberofpixels)*objectsize;
plot(sampleX,sample,'k')

% Plot style
hold on
title(['Object size = ',num2str(objectsize),' µm ; TC = ',num2str(TC*10^6),' µs ; n = ',num2str(order)])
xlabel('Distance (µm)') 
ylabel('Normalized signal (A.U.)') 
limits=[0 numberofpixels*objectsize*1.2 0 1.05]; axis(limits)

% Plot input
inputX=(0:1:(length(InputSample)-1))*objectsize/samplesperpixel;
plot(inputX,InputSample,'b')

% Plot output
plot(inputX,bestoutput{Index},'r')
legend('Object','Input','Output')
hold off

end