%Generate Fluorescently Labeled or Tagged Images of Molecules
%   (with a non-homogenous background)

%The program is used to simulate microscopic images in fluorescent 
%imaging. It can be used to verify algorithms in single particle tracking and 
%detection as well as other applications. It features a non-homogenous
%background determined through a trigonometric function that has varying
%parameters (ie. amplitude, frequency, phase shift,etc.) to generate the
%noise image. Ihe image is then looped to introduce increased
%irregularity in the noise image. Normally distributed noise is also introduced. The molecules are then appended onto this
%Molecules are then appended onto the image as per the usual method (See
%below).

%**************************************************************************

%The program uses the following variables and input parameters as a basis for 
%the subsequent image generation:

%   Molecules=number of molecules in the image

%   w=width of frame
%   h=height of frame 
%   n=the final frame based on w x h

%   mI=mean intensity
%   sI=std of intensity
%       -calculated from a user input percentage of mI

%   SNR=signal-to-noise ratio

%   mB=mean background mu
%       -calculated from the SNR and mI where, mB=(mI)/(SNR+1)
%   sB=standard deviation of the background noise
%       -calculated from a user input percentage of mB

%   Mag=the magnification of the camera 

%   lambda=the wavelength of light 

%   NA=the numerical aperture 

%   Pix_Size=the size of the pixel without magnification

%   Pix_Mag=the pixel size adjusted to the magnification

%   mPSF=mu of the point spread function

%   sPSF=std of the point spread function

%**************************************************************************
%The above parameters are used to generate the final image:

%   Image.Background --> initally generates a normally distributed background 
%   noise

%   Xi,Yi --> randomly generates the centroid positions of the molecules

%   Ii --> generates intensities from a normal distribution for each
%   molecule (intensities are capped with an if statement to never exceed 1)

%   sPSFXi,sPSFYi --> generates the standard deviation of the widths in the
%   X and Y directions independently to generate non-circular molecules

%   Image.Molecules --> generates the molecular image based on the above
%   criteria

%   Image.Output --> Image.Molecules is appended onto the initial
%   Image.Background (containing just the background at first) to generate a final 
%   Image.Output with the molecules and background noise

%**************************************************************************

%enter the number of images to generate 
Images=input('Enter number of images to generate [d/f=10]:');
if isempty(Images)
    Images=10;
end

%enter rectangular size of frame n
w=input('Enter field width (pixels) [d/f=256]:');
if isempty(w)
    w=256;
end

h=input('Enter field height (pixels) [d/f=256]:');
if isempty(h)
    h=256;
end

% generate the coordinates based on the field size n
Xpos=ones(1,w)'*(1:w);
Ypos=(1:h)'*ones(1,h);

%entering the #of molecules
Molecules=input('Enter the number of molecules [d/f=30]:');
if isempty(Molecules)
    Molecules=30;
end


%intensity and subsequent std input
mI=input('Enter a mean intensity (0-1) [d/f=0.5]:');
if isempty(mI)
    mI=0.5;
end

if mI>1
    error('The mean intensity must have a value between 0 and 1');
end

if mI<0
    error('The mean intensity must have a value between 0 and 1');
end

sI=((input('Enter intensity standard deviation (%) [d/f=10%]:')/100)*mI);
if isempty(sI)
    sI=0.1*mI;
end

%input parameters for generating the widths in the x and y directions of
%each of the molecules 
Mag=input('Enter camera magnification (X) [d/f=100x]:');%camera magnification
if isempty(Mag)
    Mag=100;
end
Pix_Size=input('Enter pixel size w/o magnification (um) [d/f=16 um]:');
if isempty(Pix_Size)
    Pix_Size=16;
end

Pix_Mag=((Pix_Size)/(Mag))*1000;%normalize pixel size to magnification

lambda=input('Enter wavelength of light (nm) [d/f=720nm]:'); %enter light wavelength
if isempty(lambda)
    lambda=720;
end

NA=input('Enter numerical aperture [d/f=1.4]:'); %enter NA
if isempty(NA)
    NA=1.4;
end

mPSF=0.5*(lambda/(2*NA))/Pix_Mag;%reduce FWHM to half for the rightbound width value 
sPSF=((input('Enter standard deviation of width (%) [d/f=10%]:')/100)*mPSF);
if isempty(sPSF)
    sPSF=0.1*mPSF;
end

%generating the noise image minimum value for the noise
SNR=input('Enter S/N Ratio [d/f=3]:');
if isempty(SNR)
    SNR=3;
end

mB=(mI)/(SNR+1); %calculate noise from a given SN ratio

sB=((input('Enter standard deviation of noise (%) [d/f=10%]:')/100)*mB);
if isempty(sB)
    sB=0.1*mB;
end

%setting up the matrix of input parameters
data=zeros(1,11); 

data(:,1)=w;
data(:,2)=h;
data(:,3)=Molecules;
data(:,4)=mI;
data(:,5)=sI;
data(:,6)=Mag;
data(:,7)=Pix_Size;
data(:,8)=lambda;
data(:,9)=NA;
data(:,10)=SNR;
data(:,11)=mB(1,1);
data(:,12)=sB;

%saving the input parameters to an ascii file 
dlmwrite('mol_input_parameters.txt','Size Molecules IntensityMU IntensitySTD Mag Pix_Size Pix_Mag Lambda NA WidthMu WidthSTD SNR NoiseMu NoiseSTD', '');
dlmwrite('mol_input_parameters.txt', data, '-append','delimiter', ' ');

%preallocate the noise matrix
Image.Background=zeros(w,h);

%generate the initital parameters for the noise image
    A=0.4*mI*rand;
    if A<0.25
        A=0.25;
    end
    B=randi([1 2],1,1)*0.01;
    C=randi([26 28],1,1);
    D=0.25*mI*rand;
    if D<0.1
        D=0.1;
    end
    E=randi([1 2],1,1)*0.01;
    F=randi([46 48],1,1);

%loop the noise image trig function to create an imperfect non-homogenous
%background
for k=1:10
%the size of the polynomial spread
    x=1:1:w;
    y=1:1:h;

    [X,Y]=meshgrid(x,y);

    %generate a slight change in the function for each iteration
    
    A=A+A*0.02*rand;
    B=B+B*0.02*rand;
    C=C+C*0.02*rand;
    D=D+D*0.02*rand;
    E=E+E*0.02*rand;
    F=F+F*0.02*rand;

    %the function for varying mu 
    Z=abs(A*cos(X.*B+C)+D*sin(Y.*E-F)); %+A*sin(B*Y);

    colormap(gray)
    surf(X,Y,Z)
    view(2)  

    Z(Z<=mB)=mB;
    %append all previous images onto the same image
    Image.Background=normrnd(Z,sB,[w,h])+Image.Background;
end

%divide the noise image by the number of loops to normalize
Image.Background=Image.Background./k;
    
    for m=1:Images 

    %introducing some variance into the noise background image(maintaining the same mean)    
    Image.Output=normrnd(Image.Background,sB);
    
    %random centroid generation
    Xi=w*rand(1,Molecules);
    Yi=h*rand(1,Molecules);
        
    %generate the normally distributed random intensities for each molecule
    Ii=abs(normrnd(mI,sI,[1,Molecules]));
    %ensure that the values for the intensity do not proceed above 1 
    for q=1:Molecules
        if Ii(q)>1
            Ii(q)=1;
        end
    end
    
    %generating the width of psf standard deviations for each molecule
    sPSFXi=abs(normrnd(mPSF,sPSF,[1,Molecules]));
    sPSFYi=abs(normrnd(mPSF,sPSF,[1,Molecules]));

    
    for p=1:Molecules
        %generate the molecular image 
        Image.Molecules=Ii(p)*exp(-((((Xpos-Xi(p)).^2)/(2*sPSFXi(p)^2))+(((Ypos-Yi(p)).^2)/(2*sPSFYi(p)^2))));
        Image.Output=Image.Output+Image.Molecules;         
    end
    
    %replace all values greater than 1 with 1 in the final image
    Image.Output(Image.Output>=1)=1;
    
    %generate matrix of generated molecule data 
    data=zeros(Molecules,6);
    data(:,1)=1:Molecules;
    data(:,2)=(Xi)';
    data(:,3)=(Yi)';
    data(:,4)=(Ii)';
    data(:,5)=(sPSFXi)'; 
    data(:,6)=(sPSFYi)';

    %saving the generated output image data as an ascii file
    dlmwrite(['mol_parameters' num2str(m) '.txt'], 'spot Xcent Ycent ints WidthSTDX WidthSTDY', 'delimiter', '');
    dlmwrite(['mol_parameters' num2str(m) '.txt'], data, '-append','delimiter', ' ');
    
    %genearate the tif image 
    imwrite(uint8(255*mat2gray(Image.Output)),['mol_image' num2str(m) '.tif'],'Compression','none');

    end
    
close all

