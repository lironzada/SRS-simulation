function output = SimulationLIA(samplechoice,modulationmethod,doublephase,numberofchannels,modulationfrequency,channeldistance,pixelsize,xpixels,ypixels,Pp,Pc,TC,pixeldwelltime,order,w01x,w01y,w02x,w02y,zR1x,zR1y,zR2x,zR2y,sample_properties,Cthermal,Cshot,Csignal,Analog,CalculatingLabel,Advanced_settings,PowerMatrix,FrequencyMatrix,Fast_Frequency,PhaseMatrix)
% Written by Bart Fokker, Vrije Universiteit Amsterdam
% Current version written on 24 February 2019
% Purpose: Simulate several multiplex SRS setups
%
%\\ Input parameters \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%
% pixelsize: pixel size in um
%
% samplechoice: choice of which sample is taken
%       0) 2D bead
%       1) 1D sample
%
% modulationmethod: choice of which modulation method is chosen
%       0) single modulation on Carrier beam
%       1) double modulation on Carrier beam
%       2) double modulation, 1st on Carrier and 2nd on Probe
%
% doublephase: choice for double-phase modulation
%       0) Regular modulation
%       1) 2 channels per frequency
%
% Analog: choice for data collection method
%       0) Digital, takes only the last sample as pixel value
%       1) Analog, takes the average of the samples as pixel value
%
% numberofchannels: number of modulation and demodulation channels
% 
% modulationfrequency: modulation frequency of channel 1 in Hz
%
% channeldistance: distance between channels in Hz
% 
% xpixels: number of xpixels
%
% ypixels: number of y pixels
%
% Pp: Probe power (mW)
%
% Pc: average Carrier power (mW)
%
% TC: Time constant (s)
%
% pixeldwelltime: pixel dwell time (s)
%
% order: filter order
%
% w01x,w01y,w02x,w02y: beam waists in X and Y direction of beam 1 and 2 (um)
%
% zR1x,zR1y,zR2x,zR2y: Rayleigh ranges in X and Y direction of beam 1 and 2 (um)
%
% sample_properties: channel, radius, xoffset and yoffset
% Radius: Radius for 2D beads or object size for 1D samples (um)
%
% xcentre: center of bead (2D) or offset of object (1D) in x direction in um
%
% ycentre: center of bead (2D) or offset of object (1D) in y direction in um
%
% Cthermal: electronic noise coefficient
%
% Cshot: shot noise coefficient
%
% Csignal: signal coefficient
%
% Analog: Digital (1 sample) or Analog (averaging) data collection
%
% For the advanced settings
% Advanced_settings: simple or advanced settings
%
% PowerMatrix: matrix containing the powers for each channel
%
% FrequencyMatrix: matrix containing the modulation frequency for each channel
%
% Fast_Frequency: The fast modulation frequency for the double modulation approach
%
% PhaseMatrix: Matrix containing the phase for each channel
%
% CalculatingLabel: Interface text showing the simulation progress
%
%\\ Output parameters \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%
% output: digital sample output
%
% Plots the simulated measurements as well
%
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% Show progress message
CalculatingLabel.Text='Calculating...';
pause(0.00000000001)

%\\ Parameter initialization and verify inputs \\\\\\\\\\\\\\\\\\\\\\\\\\\
%
% Verify input values and initialize secondary variables

% Make sure that the modulation method selection is valid
modulationmethod = round(modulationmethod);
if modulationmethod>2
    modulationmethod = 2;
elseif modulationmethod<0
    modulationmethod = 0;
end

% Make sure double phase option is valid
doublephase=round(doublephase);
if doublephase>1
    doublephase = 1;
elseif doublephase<0
    doublephase = 0;
end

% When a 1D sample is selected, the #pixels in the y direction is 1
if not(samplechoice==0)
    ypixels = 1;
end

% Initialize parameters
samplerate = 210*10^6;     % samples per second
TS = 1/samplerate;
numberofpixels = xpixels;
samplesperpixel = round(pixeldwelltime*samplerate);
samplesperum = samplesperpixel/pixelsize;

%\\ Sample creation \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%
% First, define beamshape formulas
% Secondly, define beam precision such that not every sample needs to be 
% convolved with the beam matrix. This increases the simulation speed and
% reduces the amount of memory needed.
% Thirdly, create the condensed binary sample
% Fourthly, convolve the condensed binary sample with the beam to obtain
% the condensed sample.

% Non-integrated SRS beam for 2D sample
SRSbeam = @(x,y,z,w01x,w02x,w01y,w02y,zR1x,zR2x,zR1y,zR2y)exp(-2*(x.^2./(w01x^2*(1+z.^2/zR1x^2))+y.^2./(w01y^2*(1+z.^2/zR1y^2)))-2*(x.^2./(w02x^2*(1+z.^2/zR2x^2))+y.^2./(w02y^2*(1+z.^2/zR2y^2)))).*(1./(1+z.^2/zR1x^2)+1./(1+z.^2/zR1y^2)).*(1./(1+z.^2/zR2x^2)+1./(1+z.^2/zR2y^2));

% Integrated SRS beam in the y direction for 1D samples.
% The integration was done in Mathematica.
IntySRSbeam = @(x,z,w01x,w02x,w01y,w02y,zR1x,zR2x,zR1y,zR2y)exp(-2*x.^2.*(zR1x^2./(w01x^2*(z.^2+zR1x^2))+zR2x^2./(w02x^2.*(z.^2+zR2x^2))))*sqrt(pi/2)*w01y*w02y*zR1x*zR1y*zR2x*zR2y./sqrt((z.^2+zR1x^2).*(z.^2+zR2x^2).*(w01y^2*(z.^2+zR1y^2)*zR2y^2+w02y^2*(z.^2+zR2y^2)*zR1y^2));

% Size of the beam array dependents on the SRS waist, since the array needs
% to be large enough to represent the beam while being not too large, since
% that slows the simulation down. 2 um was sufficient for the default waist.
SizeOfBeamArrayum = 2*w01x*w02x/sqrt(w01x^2+w02x^2)/0.337;

% Determine the beam precision and the x sample size of the convolution 
beamprecisionpre = 0.005*(pixelsize*(pixelsize>1)+(pixelsize<=1));
SamplesPerConvolutionStep = round(beamprecisionpre*samplesperum);
beamprecision = SamplesPerConvolutionStep/samplesperum; % um/convolution steps
SamplesPerConvolutionStepX = round(beamprecisionpre*10*samplesperum);% the x directon needs less conv samples due to interpolation
beamprecisionX = SamplesPerConvolutionStepX/samplesperum;
NumberOfConvSamples = round(pixelsize/beamprecisionX)*(numberofpixels+1);

% Create sample in binary form
if samplechoice==0
    
    % 2D bead
    
    % Sample grid
    xi2 = [fliplr(-beamprecisionX:-beamprecisionX:-SizeOfBeamArrayum/2),0:beamprecisionX:SizeOfBeamArrayum/2];
    yi2 = [fliplr(-beamprecision:-beamprecision:-SizeOfBeamArrayum/2),0:beamprecision:SizeOfBeamArrayum/2]';
    ylength = length(yi2);
    
    % Create beam matrix, assuming a homogenuous sample in the z-direction
    beammatrix = integral(@(z)SRSbeam(xi2,yi2,z,w01x,w02x,w01y,w02y,zR1x,zR2x,zR1y,zR2y),-1000,1000,'ArrayValued',1)';

    % Allocate space for the condensed binary sample (2D bead)
    NumberOfConvSamplesy = round(pixelsize/beamprecision)*ypixels+ylength;
    x = beamprecisionX:beamprecisionX:beamprecisionX*NumberOfConvSamples;
    y = (beamprecision*(1-round(ylength/2)):beamprecision:beamprecision*(NumberOfConvSamplesy-round(ylength/2)))';
    sample=false(length(y),length(x),numberofchannels);
    
    % Produce 2D bead with centre at (x0,y0)
    for i=1:numberofchannels
        properties=sample_properties{i};
        if ~isempty(properties)
            r = properties(2,1,:);
            x0 = properties(3,1,:);
            y0 = properties(4,1,:);
            sample(:,:,i) = sum(((x-x0).^2+(y-y0).^2)<=r.^2,3)>=1;
        end
    end
    
else
    
    % 1D particles
    
    % Sample x coordinates
    xi = (-SizeOfBeamArrayum/2:beamprecisionX:SizeOfBeamArrayum/2)';
    
    % Create beam matrix, assuming a homogenuous sample in the y,z-direction
    beammatrix = integral(@(z)IntySRSbeam(xi,z,w01x,w02x,w01y,w02y,zR1x,zR2x,zR1y,zR2y),-1000,1000,'ArrayValued',1)';
    
    % Allocate space for the condensed binary sample (1D sample)
    x = beamprecisionX:beamprecisionX:beamprecisionX*NumberOfConvSamples;    
    sample=false(1,length(x),1,numberofchannels);
    
    % Produce 1D sample
    for i=1:numberofchannels
        properties=sample_properties{i};
        if ~isempty(properties)
            r = properties(2,1,:);
            x0 = properties(3,1,:);
            sample(1,:,i) = sum((x>=x0)&(x<=x0+r),3)>=1;
        end
    end
end

% Sum of the beam matrix for signal normalization
beammatrixsum=sum(sum(beammatrix));

% Convolve binary sample with beam
if samplechoice==0
    
    % 2D sample
    convoutput = zeros(ypixels,NumberOfConvSamples,numberofchannels);
    
    for channeli = 1:numberofchannels
        
        % Convolve each y slice independently
        for i=1:ylength
            convoutput(:,:,channeli) = convoutput(:,:,channeli)+conv2(1,beammatrix(:,i),sample(round(i:(NumberOfConvSamplesy-ylength-1)/(ypixels-1):NumberOfConvSamplesy-ylength+i),:,channeli),'same')/beammatrixsum;
        end
        
    end
    
else
    
    % 1D sample
    convoutput = zeros(1,NumberOfConvSamples,numberofchannels);
    for channeli = 1 : numberofchannels
        convoutput(:,:,channeli) = conv(sample(:,:,channeli),beammatrix,'same')/beammatrixsum;
    end
    
end

%\\ Modulation, filtering and demodulation \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%
% First, initialize the frequencies and phases for the simple and advanced
% settings cases.
% Secondly allocate the sample output.
% Thirdly, enter a for loop which: (1) expands the condensed sample to its 
% actual size 1 pixel column (X) at a time, (2) modulates the Probe or 
% Carrier beam with the selected modulation method, (3) adds shot and
% electronic noise and adds everything together to a 2D or 1D sample array,
% (4) enters this array into the LIA simulation for demodulation and 
% filtering, (5) saves output in allocated array.

% Initialize the frequency and phase arrays for channel encoding
if Advanced_settings
    % Advanced settings: set frequencies
    Fmod=FrequencyMatrix;
    phase=PhaseMatrix;
    if modulationmethod~=0
        modulationfrequency=Fast_Frequency;
    else
        modulationfrequency=0;
    end
else
    % Simple settings
    % Set frequencies
    if doublephase==0
        Fmod=(100*10^3-(numberofchannels-1)*channeldistance:channeldistance:100*10^3)';
        phase=zeros(numberofchannels,1);
    else
        Fmod=(100*10^3-(ceil((1:numberofchannels)/2)-1)*channeldistance)';
        phase=pi/2*not(mod((1:numberofchannels),2))';
    end
    Fmod=permute(Fmod,[2,3,4,1]);
    phase=permute(phase,[2,3,4,1]);

    % Set the power per channel for the simple settings
    if modulationmethod<2
        PowerMatrix = Pc/numberofchannels*ones(1,1,1,numberofchannels);
    else
        PowerMatrix = Pp/numberofchannels*ones(1,1,1,numberofchannels);
    end
end

% Demodulation frequencies
Fdemod=unique(Fmod,'stable')+modulationfrequency;

% Demodulation frequencies
demod = permute(Fdemod,[2,3,4,1]);

% Inititialize filter output for filter start-end values and image grid
LastSample = zeros(ypixels,1,order+1,length(Fdemod));
DemodulatedSample = zeros(ypixels,numberofpixels,1,length(Fdemod));

% Initialize sample parameters for keeping track of the sample size
ConvSamplesPerPixel = round(pixelsize/beamprecisionX);
SamplesPerPixel = ConvSamplesPerPixel*SamplesPerConvolutionStepX;

% Initialize the demodulation array
demodnoise = zeros(ypixels,SamplesPerPixel+1,order+1,length(Fdemod));
    
% for-loop for demodulating each sequential pixel in the x direction. The 
% for-loop is used because it prevents the matrix sizes from becoming too
% large.
for pixel=1:numberofpixels
    
    % Expanding the sample array to its true size
    InputSample=zeros(ypixels,SamplesPerPixel,numberofchannels);
    if samplechoice==0 && ypixels>1
        for i=1:numberofchannels
            InputSample(:,:,i)=interp2(convoutput(:,1+(pixel-1)*ConvSamplesPerPixel:pixel*ConvSamplesPerPixel+1,i),1:1/SamplesPerConvolutionStepX:ConvSamplesPerPixel+1-1/SamplesPerConvolutionStepX,(1:size(convoutput,1))');
        end
    else
        for i=1:numberofchannels
            InputSample(:,:,i)=interp1(convoutput(1,1+(pixel-1)*ConvSamplesPerPixel:pixel*ConvSamplesPerPixel+1,i),1:1/SamplesPerConvolutionStepX:ConvSamplesPerPixel+1-1/SamplesPerConvolutionStepX);
        end
    end
    
    % Add material channels together.
    InputSample2=sum(InputSample.*permute(Csignal,[3,4,2,1]),3);
    
    % Indices for all samples in this for-loop entry
    SampleIndices = 1+(pixel-1)*SamplesPerPixel:pixel*SamplesPerPixel;
    
    % Determine Carrier and Probe power depending on the modulation method
    if modulationmethod==0
        
        % Split Carrier times oscillations
        ProbePower = Pp;
        CarrierPower = sum(PowerMatrix.*InputSample2.*(cos(2*pi*(modulationfrequency+Fmod).*(SampleIndices)*TS+phase)),4);
        
    elseif modulationmethod==1
        
        % Split Carrier times high and low frequency oscillations
        ProbePower = Pp;
        CarrierPower = 2*sum(PowerMatrix.*InputSample2.*(cos(2*pi*modulationfrequency*(SampleIndices)*TS)).*(1+cos(2*pi*Fmod.*(SampleIndices)*TS+phase)),4);
        
    else
        
        % Split Probe times low frequency oscillation and Carrier times high frequency oscillation
        ProbePower = 2*sum(PowerMatrix.*(1+cos(2*pi*Fmod.*(SampleIndices)*TS+phase)),4);
        PumpPowerG = 2*sum(PowerMatrix.*InputSample2.*(1+cos(2*pi*Fmod.*(SampleIndices)*TS+phase)),4);
        CarrierPower = Pc*(cos(2*pi*modulationfrequency*(SampleIndices)*TS));
        
    end

    % Determine the noise with sqrt(samplerate)/2 is a conversion factor.
    ShotNoise = sqrt(samplerate)/2*Cshot*sqrt(ProbePower).*wgn(ypixels,SamplesPerPixel,0);
    ElectronicNoise = sqrt(samplerate)/2*Cthermal*wgn(ypixels,SamplesPerPixel,0);
    
    % Determine the SRL (L.Z. removed factor sqrt(2))
    if modulationmethod<2
        PumpPowerGain = CarrierPower.*ProbePower;
    else
        PumpPowerGain = CarrierPower.*PumpPowerG;
    end
    
    % Input sample
    input = ShotNoise+ElectronicNoise+PumpPowerGain;
    
    % Demodulation and filtering by digital LIA, which demodulates and 
    % filters the input sample with a digital RC filter as for the HF2LI of
    % Zurich Instruments.

    % Last demodulation output of previous pixel
    demodnoise(:,1,:,:) = LastSample;

    % Demodulate input sample
    demodnoise(:,2:end,1,:) = input*sqrt(2).*exp(-1i*2*pi*demod.*SampleIndices*TS);

    % Digital RC filter
    for i=2:size(input,2)+1

        demodnoise(:,i,2:end,:)=exp(-TS./TC).*demodnoise(:,i-1,2:end,:)+(1-exp(-TS./TC)).*demodnoise(:,i-1,1:end-1,:);

    end

    % Select the last sample for each filter output of this pixel as a starting
    % value for the next pixel.
    LastSample = demodnoise(:,end,:,:);

    % Digital or Analog data collection 
    if Analog==0
        % Digital: take last sample
        sampleoutput = LastSample(:,end,end,:);
    else
        % Analog: take average of all samples
        sampleoutput = mean(demodnoise(:,:,end,:),2);
    end
    
    % Add pixel value
    DemodulatedSample(:,pixel,1,:) = sampleoutput;
    
    % Update progress message
    CalculatingLabel.Text=[num2str(round(100*pixel/numberofpixels)),'%'];
    pause(0.00000000001)
end

%\\ Create image \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%
% First, create image matrix for a gray scale image
% Secondly, plot the image of all channels in the same figure using subplots

% Convert output matrix (x,y,1,channel) to (x,y,channel) 
output = permute(DemodulatedSample,[1,2,4,3]);

% Produce grayscale image matrix
if doublephase==0
    ImageMatrix=(real(output)-min(min(real(output))))./(max(max(real(output)))-min(min(real(output))))*255;
    
else
    ImageMatrix = zeros(xpixels,ypixels,numberofchannels);
    F=unique(Fmod,'stable');
    
    for i=1:numberofchannels
        Index=find(F==Fmod(i));
        if phase(i)
            % phase = pi
            ImageMatrix(:,:,Index)=(imag(output)-min(min(imag(output))))./(max(max(imag(output)))-min(min(imag(output))))*255;
        else
            % phase = 0
            ImageMatrix(:,:,Index)=(real(output)-min(min(real(output))))./(max(max(real(output)))-min(min(real(output))))*255;
        end
    end
end

% Plot image matrix in the same figure
figure
for i=1:numberofchannels
    
    % Specify subplot
    subplot(ceil(numberofchannels/2),ceil(numberofchannels/ceil(numberofchannels/2)),i)
    
    % Plot image
    if samplechoice==0
        imshow(ImageMatrix(:,:,i),[0,255])
    else
        % Make 1D image clearer by increasing its vertical size
        imshow(repmat(ImageMatrix(:,:,i),round(xpixels/10),1),[0,255])
    end
    
end

end
