%Generate Single Frame Images of Fluorescently Emitting Molecules in a
%Microscope

%The program is used to simulate microscopic images in fluorescent 
%imaging. It can be used to verify algorithms in single particle tracking and 
%detection as well as other applications.

%**************************************************************************

%The program uses the follwing variables and input parameters as a basis for 
%the subsequent image generation:

%   Images=the number of images to generate under these parameters

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

%   Image.Output --> initally generates a normally distributed background 
%   noise

%   Xi,Yi --> randomly generates the centroid positions of the molecules

%   Ii --> generates intensities from a normal distribution for each
%   molecule (intensities are capped with an if statement to never exceed 1)

%   sPSFXi,sPSFYi --> generates the standard deviation of the widths in the
%   X and Y directions independently to generate non-circular molecules

%   Image.Molecules --> generates the molecular image based on the above
%   criteria

%   Image.Output --> Image.Molecules is appended onto the initial
%   Image.Output (containing just the background at first) to generate a final 
%   Image.Output with the molecules and background noise

%**************************************************************************

%enter number of images to generate
Images=input('Enter number of images to generate:');

%enter rectangular size of frame n
w=input('Enter field width (pixels):');
h=input('Enter field height (pixels):');

% generate the coordinates based on the field size n
X=ones(1,w)'*(1:w);
Y=(1:h)'*ones(1,h);

%entering the #of molecules
Molecules=input('Enter the number of molecules:');

%intensity and subsequent std input
mI=input('Enter a mean intensity (0-1):');

if mI>1
    error('The mean intensity must have a value between 0 and 1');
end

if mI<0
    error('The mean intensity must have a value between 0 and 1');
end

sI=((input('Enter intensity standard deviation (%):')/100)*mI);

%input parameters for generating the widths in the x and y directions of
%each of the molecules 
Mag=input('Enter magnification:');%camera magnification
Pix_Size=input('Enter pixel size w/o magnification (um):');

Pix_Mag=((Pix_Size)/(Mag))*1000;%normalize pixel size to magnification

lambda=input('Enter wavelength of light (nm):'); %enter light wavelength
NA=input('Enter numerical aperture:'); %enter NA

mPSF=0.5*(lambda/(2*NA))/Pix_Mag;
sPSF=((input('Enter standard deviation of width (%):')/100)*mPSF);

%generating the noise image
SNR=input('Enter S/N Ratio:');
mB=(mI)/(SNR+1); %calculate noise from a given SN ratio

sB=((input('Enter standard deviation of noise (%):')/100)*mB);

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
data(:,11)=mB;
data(:,12)=sB;

%saving the input parameters to an ascii file 
dlmwrite('mol_input_parameters.txt','Size Molecules IntensityMU IntensitySTD Mag Pix_Size Pix_Mag Lambda NA WidthMu WidthSTD SNR NoiseMu NoiseSTD', '');
dlmwrite('mol_input_parameters.txt', data, '-append','delimiter', ' ');

%program is looped here to create multiple images under the same parameters
for m=1:Images 
    %generating the noise image 
    Image.Output=normrnd(mB,sB,[w,h]);

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
        Image.Molecules=Ii(p)*exp(-((((X-Xi(p)).^2)/(2*sPSFXi(p)^2))+(((Y-Yi(p)).^2)/(2*sPSFYi(p)^2))));
        Image.Output=Image.Output+Image.Molecules; 
        %genearate the tif image 
        imwrite(uint8(255*mat2gray(Image.Output)),['mol_image' num2str(m) '.tif'],'Compression','none');

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
        
    end

end
