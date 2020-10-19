function [OptimalTC,CNRproportional,responseoutput,bestoutput,sample,InputSample,samplesperpixel,numberofpixels,beamprecision]=LIAOptimalSettings(objectsize,choice,order,TC,w01x,w01y,w02x,w02y,zR1x,zR2x,zR1y,zR2y)
% Written by Bart Fokker, Vrije Universiteit Amsterdam
% written:10 January 2019
% Edited by Liron Zada
% Last Edit: 14 May 2020
% Purpose: Find the settings for which the CNR is maximal.
% The sample is a serie of block functions or a single block. The optimal
% settings were determined by checking for which TC the filter responds 
% best; results in the largest CNR. In this code the objectsize and filter
% order are the free parameters. The input TCs are only meant to find the
% optimal TC for each specific setting. The object size functions here as 
% the resolution and the pixel size. The object size is the length of the 
% objects and the distance between them in um.
%
%\\ Input parameters \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%
% objectsize: object size in um
%
% choice: single particle or tightly packed particles
%       0) Tightly packed particles
%       1) Single particle
%
% order: filter orders
%
% TC: TC range
% 
% w01x,w01y,w02x,w02y,zR1x,zR2x,zR1y,zR2y: beam shape properties
%
%\\ Output parameters \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%
% OptimalTC: Optimal TC
% 
% CNRproportional: CNR proportional value for the optimal TC
%
% responseoutput: Peak (- Valley) for the optimal TC
%
% Plotting outputs:
%
%   - bestoutput: Demodulation output for the optimal filter setting
%
%   - sample: Object (top hat or block wave)
%
%   - InputSample: Input signal (convolution sample with SRS focal volume)
%
%   - samplesperpixel: Samples per pixel
%
%   - numberofpixels: Number of pixels
%
%   - beamprecision: precision of beam in um (for the convolution)
%
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% Initialize parameters
samplerate=210*10^6; % samples per second
TS=1/samplerate;
pixeldwelltime=100.85*10^-6; % in seconds
numberofpixels=10;
samplesperpixel=round(pixeldwelltime*samplerate);
samplesperum=samplesperpixel/objectsize;

%\\ Sample creation \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%
% First, define beamshape formulas
% Secondly, define beam precision such that not every sample needs to be 
% convolved with the beam matrix. This increases the simulation speed and
% reduces the amount of memory needed.
% Thirdly, create the condensed pre-sample
% Fourthly, convolve the condensed pre-sample with the beam
% Fifthly, expand pre-sample into input sample
%
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% The y direction is already integrated.
% Integrate in the x direction.
IntySRSbeam=@(x,z,w01x,w02x,w01y,w02y,zR1x,zR2x,zR1y,zR2y)exp(-2*x.^2.*(zR1x^2./(w01x^2*(z.^2+zR1x^2))+zR2x^2./(w02x^2.*(z.^2+zR2x^2))))*sqrt(pi/2)*w01y*w02y*zR1x*zR1y*zR2x*zR2y./sqrt((z.^2+zR1x^2).*(z.^2+zR2x^2).*(w01y^2*(z.^2+zR1y^2)*zR2y^2+w02y^2*(z.^2+zR2y^2)*zR1y^2));

% Size of the beam array is chosen as 10 um.
SizeOfBeamArrayum=10;

% Determine the beam precision and the x sample size of the convolution 
SamplesPerConvolutionStep=round(0.005*samplesperum);
beamprecision=SamplesPerConvolutionStep/samplesperum;
xi=(-SizeOfBeamArrayum/2:beamprecision:SizeOfBeamArrayum/2)';
beamarray=integral(@(z)IntySRSbeam(xi,z,w01x,w02x,w01y,w02y,zR1x,zR2x,zR1y,zR2y),-1000,1000,'ArrayValued',1)';

% Pre-sample building blocks
sample0=zeros(1,round(objectsize/beamprecision));
sample1=ones(1,round(objectsize/beamprecision));

% Create pre-sample
if choice==1
    % Tightly packed particles
    sample=[sample0,sample0,repmat([sample1,sample0],1,round((numberofpixels-2)/2))];
    
else
    % Single particle
    sample=[sample0,sample0,sample1,repmat(sample0,1,numberofpixels-3)];
    
end

% Normalize and convolve sample
samplex=sample/sum(beamarray);
convoutput=conv(samplex,beamarray,'same');

% Expand sample into input values
InputSample=interp1(convoutput,1:1/SamplesPerConvolutionStep:length(samplex)+1-1/SamplesPerConvolutionStep);

%\\ Filter \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%
% Firstly: Obtain response through the digital LIA simulation
% Secondly: Obtain through this the maximal CNR
% Thirdly: Select through this the optimal TC
% Fourthly: Repeat previous steps with a greater TC precision
%
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% Filter the input sample
index=-1;
while (index==1||index==length(TC)||index==-1)
    
    % If the optimal TC is the last or first entry, it likely means that
    % the optimal TC lies outside its range. Therefore, run the simulation
    % again with a larger range.
    if (index==1)
        TC=(max([TC(index)-20*10^-6,10^-6]):10^-6:TC(index)+10^-6)';
    elseif (index==length(TC))
        TC=(TC(index)-10^-6:10^-6:TC(index)+20*10^-6)';
    end
    
    % Obtain response through LIA simulation
    response1=TCselection(TC,order,length(InputSample)/numberofpixels,numberofpixels,TS,InputSample,choice);

    % Calculate niose from the filter order and the TC 
    NEPBW=gamma(order-1/2)./(4*sqrt(pi).*TC.*gamma(order));
    NoiseProportional1=1.*sqrt(NEPBW);% This is not the true noise.
    
    % Obtain the CNR and the TC for which the CNR is maximal
    CNRProportional1=response1./NoiseProportional1;
    [~,index]=max(CNRProportional1);
    
    % In case the optimal TC is smaller than 1 us, stop the loop to prevent
    % an infinite loop
    if (length(index)==2)
        break;
    end
end

% Repeat simulation with greater precission
TC2=TC(index)+10^-6*(-0.9:0.05:0.9)';
[response2,output]=TCselection(TC2,order,length(InputSample)/numberofpixels,numberofpixels,TS,InputSample,choice);

% Calculate from order and TC the noise.
NEPBW=gamma(order-1/2)./(4*sqrt(pi).*TC2.*gamma(order));
NoiseProportional2=1.*sqrt(NEPBW);% This is not the true noise.

% Obtain optimal TC with its response and proportional CNR value
CNRProportional2=response2./NoiseProportional2;
[CNRproportional,index2]=max(CNRProportional2);
OptimalTC=TC2(index2);
responseoutput=response2(index2);

% Select output of best filter settings
bestoutput=output(index2,:);

end


function [response,output]=TCselection(TC,k,samplesperpixel,numberofpixels,TS,input,choice)
% Digital LIA simulation function
%
% Input
% TC:               Time constants
% k:                filter order
% samplesperpixel:  samples per pixel
% numberofpixels:   number of pixels
% TS:               sample time
% input:            input sample
% choice:           single or tightly packed particles
%
% Output
% response:         maximal response or contrast
% output:           demodulation output for plot
%
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% Sample length
llength=length(TC);

% Initialize the demodulation array
demodnoise=zeros(llength,round(samplesperpixel*numberofpixels),k+1);

% Demodulate input sample
demodnoise(:,1:end,1,:) = input.*ones(llength,1);

% Digital RC filter
for i=2:size(input,2)+1

    demodnoise(:,i,2:end,:)=exp(-TS./TC).*demodnoise(:,i-1,2:end,:)+(1-exp(-TS./TC)).*demodnoise(:,i-1,1:end-1,:);

end

% Select the last sample for each filter output to use in the plot
LIAoutput = demodnoise(:,:,k+1,:);

% Output for plot
output=real(LIAoutput(:,2:1:end));

% Obtain response for optimal setting
if choice==0

    % For single particle: response is maximum value
    response=max(LIAoutput,[],2);

else

    % For tightly packed particles: response is contrast between maximum
    % and minimum value.

    % Discrete derivative 
    dx=LIAoutput(:,2:end)-LIAoutput(:,1:end-1);

    % Obtain maxima and minima
    positivematrix=dx(:,2:end)<0&dx(:,1:end-1)>0;
    negativematrix=dx(:,1:end-1)<0&dx(:,2:end)>0;

    % Find first maximum and minimum
    [~,positiveindex]=max(positivematrix,[],2);
    [~,negativeindex]=max(negativematrix,[],2);
    LinearPositiveIndices=sub2ind(size(LIAoutput),(1:llength)',positiveindex+1);
    LinearNegativeIndices=sub2ind(size(LIAoutput),(1:llength)',negativeindex+1);

    % Delete first maximum and minimum
    positivematrix(LinearPositiveIndices-llength)=0;
    negativematrix(LinearNegativeIndices-llength)=0;

    % Find second maximum and minimum
    [~,positiveindex]=max(positivematrix,[],2);
    [~,negativeindex]=max(negativematrix,[],2);
    LinearPositiveIndices=sub2ind(size(LIAoutput),(1:llength)',positiveindex+1);
    LinearNegativeIndices=sub2ind(size(LIAoutput),(1:llength)',negativeindex+1);

    % Obtain second maximum and minimum
    positiveresponse=LIAoutput(LinearPositiveIndices);
    negativeresponse=LIAoutput(LinearNegativeIndices);

    % Obtain response (here as contrast)
    responsecheck=(samplesperpixel*2>LinearPositiveIndices-LinearNegativeIndices)&(LinearNegativeIndices-LinearPositiveIndices>0);
    response=(positiveresponse-negativeresponse).*responsecheck;

end

end
