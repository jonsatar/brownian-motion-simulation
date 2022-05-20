%Particle Simulation of Brownian Motion with Trajectories

%The program is used to simulate brownian motion in microscopic images 
%in fluorescent imaging. It can be used to verify algorithms in single
%particle tracking as well as other applications. A series of still images can be
%generated for molecules within a frame, its corresponding data, and
%concatenation into a movie file. The trajectories of the particles are
%also generated for these purposes.

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

%   t=time interval between frames

%   T=total time interval 

%   A diffusion coefficient is entered and then vectorized to guide the movement of the
%   molecules in the frame. These changes are then appended onto the
%   original displacements and looped to create a certain number of frames
%   generating this movement multiple times in subsequent images. Also, the 
%   standard deviation of each
%   individual molecule is halved from the user input percentage and there
%   is allowed to be some variation in the intensities, the sigmaX and sigmaY
%   upon each iteration or molecule movement. 
%**************************************************************************

%The above parameters are used to generate the final image described below:

%   Image.Output --> initally generates a normally distributed background 
%   noise

%   Xi,Yi --> randomly generates the centroid positions of the molecules

%   Ii --> generates intensities from a normal distribution for each
%   molecule (intensities are capped with an if statement to never exceed 1
%   as well)

%   sPSFXi,sPSFYi --> generates the standard deviation of the widths in the
%   X and Y directions independently to generate non-circular molecules

%   Image.Molecules --> generates the molecular image based on the above
%   criteria

%   Image.Output --> Image.Molecules is appended onto the initial
%   Image.Output (containing just the background at first) to generate a final 
%   Image.Output with the molecules and background noise

%**************************************************************************

%enter size
w=input('Enter width of the frame (pixels):');
h=input('Enter height of the frame (pixels):');

% generate the coordinate matrix
X=ones(1,w)'*(1:w);
Y=(1:h)'*ones(1,h);

%#of molecules
Molecules=input('Enter number of molecules:');

%intensity and subsequent std input
mI=input('Enter mean intensity (0-1):');

%error condition to ensure a correct intensity value between 0 and 1
if (mI>1)
    error('The mean intensity must have a value between 0 and 1');
end
if (mI<0)
    error('The mean intensity must have a value between 0 and 1');
end

%calculate a standard deviation of the intensity
sI=((input('Enter standard deviation of the intensity (%):')/100)*mI);

%input parameters for generating the widths in the x and y directions of
%each of the molecules 
Mag=input('Enter camera magnification (X):');%camera magnification
Pix_Size=input('Enter pixel size without magnification (um):');
Pix_Mag=((Pix_Size)/(Mag))*1000;%normalize pixel size to magnification
Lambda=input('Enter wavelength of light (nm):'); %enter light wavelength
NA=input('Enter numerical aperture:'); %enter NA

mPSF=0.5*((Lambda/(2*NA))/Pix_Mag);%multiply the FWHM by 0.5
sPSF=(input('Enter standard deviation of width (%):')/100)*mPSF; %std of width calculation

%generating the criteria for the noise image construction
SNR=input('Enter S/N ratio:');
mB=(mI)/(SNR+1); %calculate noise from a given SN ratio
sB=(input('Enter standard deviation of the noise (%):')/100)*mB;

%generate the noise image based off of the above parameters
Image.Output=normrnd(mB,sB,[w,h]);

%setting up the matrix of input parameters
data=zeros(1,11); 
data(:,1)=w;
data(:,2)=h;
data(:,3)=Molecules;
data(:,4)=mI;
data(:,5)=sI;
data(:,6)=Mag;
data(:,7)=Pix_Size;
data(:,8)=Lambda;
data(:,9)=NA;
data(:,10)=SNR;
data(:,11)=mB;
data(:,12)=sB;

%saving the input parameters to an ascii file 
dlmwrite('mol_input_parameters.txt','width height Molecules IntensityMU IntensitySTD Mag Pix_Size Pix_Mag Lambda NA SNR NoiseMu NoiseSTD', '');
dlmwrite('mol_input_parameters.txt', data, '-append','delimiter', ' ');

% generate random centroids for each molecule
Xi=w*rand(1,Molecules);
Yi=h*rand(1,Molecules);

%generate random intensities for each molecule from a gaussian PSF
Ii=abs(normrnd(mI,sI,[1,Molecules]));

%if the intensity of the molecule achieves a value greater than 1, then
%reduce its value to 1 
for q=1:Molecules
if Ii(q)>1
Ii(q)=1;
end
end

%generate the width of the psf standard deviations
sPSFXi=abs(normrnd(mPSF,sPSF,[1,Molecules]));
sPSFYi=abs(normrnd(mPSF,sPSF,[1,Molecules]));

%generate the molecular output image with the noise background
for p=1:Molecules
Image.Molecules=Ii(p)*exp(-((((X-Xi(p)).^2)/(2*sPSFXi(p)^2))+(((Y-Yi(p)).^2)/(2*sPSFYi(p)^2))));
Image.Output=Image.Output+Image.Molecules;
end

%save a tif image 
imwrite(uint8(255*mat2gray(Image.Output)),'mol_imageORIG.tif','Compression','none');

%generate matrix of generated molecule data 
data=zeros(Molecules,6);
data(:,1)=1:Molecules;
data(:,2)=(Xi)';
data(:,3)=(Yi)';
data(:,4)=(Ii)';
data(:,5)=(sPSFXi)'; 
data(:,6)=(sPSFYi)';

%save a variable of the original molecular data to append the later dI, dSx
%and dSy values onto for consistent variation
IOriginal=(Ii)';
SxOriginal=(sPSFXi)';
SyOriginal=(sPSFYi)';

%saving the generated output image data as an ascii file
dlmwrite('mol_parametersORIG.txt', 'spot Xcent Ycent ints WidthSTDX WidthSTDY', 'delimiter', '');
dlmwrite('mol_parametersORIG.txt', data, '-append','delimiter', ' ');

%enter the criteria for creating looped multiple random walks 
t=input('Enter time between frames (s):');
T=input('Enter desired time interval (s):');
Num_Frames=round(T/t);
D_um=input('Enter diffusion coefficient (um^2/s):');
D_pix=((D_um)*1000^2)/(Mag^2); %normalize diffusion coefficient to pix^2/s

MSD=4*D_pix*t; %calculate the mean square displacement

%preallocate matrix for speed, and generate the particle and frame entries
%in the track data matrix
trackdata=zeros(Num_Frames*Molecules,7);
for n=0:(Molecules-1)
trackdata((1+n*Num_Frames):((n+1)*Num_Frames),1)=1:Num_Frames;
trackdata((1+n*Num_Frames):((n+1)*Num_Frames),2)=repmat(n+1,1,Num_Frames);
end

%generating a new image consistent with the MSD parameters
for m=1:Num_Frames
 
%generate a new noise background image for each iteration
Image.Output=normrnd(mB,sB,[w,h]);

    for p=1:Molecules
    
    %generate a radius of displacement based off of a normrnd distribution
    r=normrnd(0,sqrt(MSD));  
    %generate a random angle for the displacement (multiply by pi only to
    %generate an ange in radians that is between 0-180 in degrees)
    theta=rand*pi;
    %calculate the x and y displcements
    dx=r*cos(theta);
    dy=r*sin(theta);
    %append the change of x and y onto the previous positions
    data(p,2)=data(p,2)+dx;
    data(p,3)=data(p,3)+dy;
    
    %vary the itensity for each molecule slightly
    dI=normrnd(0,0.5*sI);%note reduced std
    data(p,4)=IOriginal(p,1)+dI;
    
    %if the intensity achieves a value greater than 1, reduce its value to
    %1 
    if data(p,4)>1
        data(p,4)=1;
    end
    
    %generate standard deviations that fluctuate when the molecules move 
    dSx=normrnd(0,0.5*sPSF);
    data(p,5)=SxOriginal(p,1)+dSx;
    
    dSy=normrnd(0,0.5*sPSF);
    data(p,6)=SyOriginal(p,1)+dSy;
    
    %calculate a new image based off the above changes 
    Image.Molecules=data(p,4)*exp(-((((X-data(p,2)).^2)/(2*data(p,5)^2))+(((Y-data(p,3)).^2)/(2*data(p,6)^2))));
    Image.Output=Image.Output+Image.Molecules;   
 
    end
    
    %saving parallel molecular data into a trackdata matrix
    trackdata((1+(m-1)):Num_Frames:Num_Frames*Molecules,3)=data(:,2); %saving the x data
    trackdata((1+(m-1)):Num_Frames:Num_Frames*Molecules,4)=data(:,3); %saving the y data 
    trackdata((1+(m-1)):Num_Frames:Num_Frames*Molecules,5)=data(:,4);%saving the intensity
    trackdata((1+(m-1)):Num_Frames:Num_Frames*Molecules,6)=data(:,5);%saving the sigmaX
    trackdata((1+(m-1)):Num_Frames:Num_Frames*Molecules,7)=data(:,6);%saving the sigmaY
    
    dlmwrite('mol_track_parameters.txt', 'Frame Molecule Xcent Ycent Intensity WidthSTDX WidthSTDY', 'delimiter', '');
    dlmwrite('mol_track_parameters.txt', trackdata, '-append','delimiter', ' ');
    
    %save the new images data to an ascii file
dlmwrite(['mol_parameters' num2str(m) '.txt'], 'Molecule Xcent Ycent Intensity WidthSTDX WidthSTDY', 'delimiter', '');
dlmwrite(['mol_parameters' num2str(m) '.txt'], data, '-append','delimiter', ' ');
  %save the new image as a tif file
imwrite(uint8(255*mat2gray(Image.Output)),['mol_image' num2str(m) '.tif'],'Compression','none'); %save image as a tiff file (8 bit output image)   

%generate the grayscale movie image frames
mov(m)=im2frame(repmat(uint8(255*mat2gray(Image.Output)), [1 1 3])); %#ok<SAGROW>

end

fps=input('Enter the number of frames per second:');

%save the movie frames to an avi file at 1 frame per second
movie2avi(mov, 'mol_movie','Compression','none', 'fps',fps);


%generate the trajectory files and then save them as .jpg
for h=1:Molecules
molecule_traj=plot(trackdata(((h-1)*Num_Frames+1):(h*Num_Frames),3),trackdata(((h-1)*Num_Frames+1):(h*Num_Frames),4),'-o');%#ok<MSNU,MSNU>
 xlabel('X-Location')
 ylabel('Y-Location')
 title(['Molecular Trajectory, Molecule: ' num2str(trackdata(((h-1)*Num_Frames+1),2))])
saveas(molecule_traj,['molecule_traj' num2str(h) '.jpg'])
end