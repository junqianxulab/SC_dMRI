function ComputeWriteDDEcordImages(niftifilein, Bvalmax, nB0, nReps)
% function ComputeWriteDDEcordImages(fname, Bvalmax, nB0, nReps)
%
% Function to create maps of axial diffusivity from filtered Diffusion
% (double diffusion type experiments) specific to a spinal cord protocol
%
% Inputs:
%       niftifilein - filename of nifti image containing diffusion weighted
%                   images.  The 4th dimension must equal the b-table below (20 dirs)
%
%       Bvalmax - Bvalue maximum set at the scanner. This is multiplied by
%       the square of the vectors to get the b-value along each direction
%       (parallel and perpendicular. Required, no default.
%
%       nReps - number of repetitions of the entire set (b0 + dirs)
%               Defaults to 1 if not included.
%
%       nB0 - number of non-diffusion weighted images (b=0).  
%               Defaults to 1 if not included.
%
% Outputs:
%   Nifti files are written to the same directory with a naming scheme
%   The prefix is same as input, but with _fDWI_Metic.nii as output.
%   Primarily, Daxial and Dradial are quantitative diffusivities
%   MeanRad is the mean high b-value (perpendicular) image.
%
%  Uses the Nifti file reader/writer functions from Matlab central:
%  https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image?requestedDomain=www.mathworks.com
%
%  Matthew Budde, PhD, Medical College of Wisconsin, 10/2020
%

    % Siemens Vector fomatting.
    % Vector [0] = (0,0.8160,0)
    % Vector [1] = (0,0.8160,0.2890)
    % Vector [2] = (0,0.8160,0.4080)
    % Vector [3] = (0,0.8160,0.5000)
    % Vector [4] = (0,0.8160,0.5770)
    % Vector [5] = (0,-0.8160,0)
    % Vector [6] = (0,-0.8160,0.2890)
    % Vector [7] = (0,-0.8160,0.4080)
    % Vector [8] = (0,-0.8160,0.5000)
    % Vector [9] = (0,-0.8160,0.5770)
    % Vector [10] = (0,0.8160,0)
    % Vector [11] = (0,0.8160,-0.2890)
    % Vector [12] = (0,0.8160,-0.4080)
    % Vector [13] = (0,0.8160,-0.5000)
    % Vector [14] = (0,0.8160,-0.5770)
    % Vector [15] = (0,-0.8160,0)
    % Vector [16] = (0,-0.8160,-0.2890)
    % Vector [17] = (0,-0.8160,-0.4080)
    % Vector [18] = (0,-0.8160,-0.5000)
    % Vector [19] = (0,-0.8160,-0.5770)

   %these values are hard-coded here for a specific protocol 
    Vector = [
    [0,0.8160,0];
    [0,0.8160,0.2890];
    [0,0.8160,0.4080];
    [0,0.8160,0.5000];
    [0,0.8160,0.5770];
    [0,-0.8160,0];
    [0,-0.8160,0.2890];
    [0,-0.8160,0.4080];
    [0,-0.8160,0.5000];
    [0,-0.8160,0.5770];
    [0,0.8160,0];
    [0,0.8160,-0.2890];
    [0,0.8160,-0.4080];
    [0,0.8160,-0.5000];
    [0,0.8160,-0.5770];
    [0,-0.8160,0];
    [0,-0.8160,-0.2890];
    [0,-0.8160,-0.4080];
    [0,-0.8160,-0.5000];
    [0,-0.8160,-0.5770]
    ];

if ~exist('Bvalmax','var')
    error('Requires max bvalue as second argument');
end

if ~exist('nB0','var')
    nB0 = 1;
    disp('Assuming 1 b=0 image; include 3rd argument if incorrect');
end

if ~exist('nReps','var')
    nReps = 1;
    disp('Assuming 1 repetition; include 4th argument if incorrect');
end



    % flags for the calculations of Daxial, Dradial.
    %Daxial uses only the AxOn (filter) images (and b-values from 3rd column
    %Dradial uses the perpendicular only
    %The filter-only images are output as MeanRad (3rd column is zero, 2nd column is max

            FiltFinal = [zeros(1,nB0) 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0];
            FiltOnly =  [ones(1,nB0)  1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0];
            AxOn =      [zeros(1,nB0) 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
            B0Ind =     [ones(1,nB0)  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
            
            AxBvals = [zeros(1,nB0) Vector(:,3)'].^2 .* Bvalmax;
            RadBvals = [zeros(1,nB0) Vector(:,2)'].^2 .* Bvalmax;

            BFilt = repmat(RadBvals,[1 nReps]);
            BAx = repmat(AxBvals,[1 nReps]);
            FiltFinal = repmat(FiltFinal,[1 nReps]);
            FiltOnly = repmat(FiltOnly,[1 nReps]);
            AxOn = repmat(AxOn,[1 nReps]);

    %setup some universal indices for the fitting.
    FiltFinalInd = find(FiltFinal);
    RadInd = find(FiltOnly);
    AxInd = find(AxOn);
    B0Ind = find(B0Ind);


    %divide bvalues by 1000 for final units of um^2/ms
    % include column of 1's for the least squares fitting to also fit the
    % offset, even though the data gets normalized
    fitAxBvals = [ones(1,size(AxInd,2)); BAx(AxInd)./-1000]';
    fitRadBvals = [ones(1,size(RadInd,2)); BFilt(RadInd)./-1000]';


    % load in the file and capture the naming prefix 
    [fp, fn, fe] = fileparts(niftifilein);
    [fp, fn, fe] = fileparts(fn);
    fe='.nii';
    OutputPrefix = strcat(fn,'_fDWI');


    % Load the data
    dwinii = load_untouch_nii(niftifilein);  %load untounch to preserve all other image header info, especially the sform/qform
    dwinii.img = double(dwinii.img); %convert to double
    dwisize = size(dwinii.img);
    
    if dwisize(4) ~= length(AxOn)
        warning('Number of images in 4th dimension do not match number of bvalues used.\n Images may process but are unlikely to be correct!');
    end

    %allocate matrices within a single structure
    Out.Daxial = zeros(dwisize(1:3));
    Out.MeanFilt = mean(dwinii.img(:,:,:,FiltFinalInd),4);
    Out.Dradial = zeros(dwisize(1:3));
    
    %dont include the S0 images in output for now
 %   Out.AxSzero = zeros(dwisize(1:3));
 %   Out.RadSzero = zeros(dwisize(1:3));

    % do the fitting.
    disp('Computing Maps...');
    for zz=1:size(dwinii.img,3)
        fprintf('Slice %d of %d',zz,size(dwinii.img,3));
        for xx = 1:size(dwinii.img,1)
            % a single dot for every 10 image lines.
            if mod(xx,10)==0
                fprintf('.');
            end
            for yy = 1:size(dwinii.img,2)

                % axial fitting using single parallel directions (with filter) and
                % monoexponential weighted least squares
                voxdat = squeeze(dwinii.img(xx,yy,zz,AxInd));
                normdat = squeeze(dwinii.img(xx,yy,zz,FiltFinalInd));
                meannormdat = mean(normdat);
                if meannormdat ~= 0
                    %normalize and log the data
                    normvdat = log(voxdat./mean(normdat));
                    %least squares
                    lls = fitAxBvals\normvdat;

                    % weighted fitting.
                    if ~any(isnan(lls))
                        w = diag(exp(fitAxBvals*lls));
                        lls = (w*fitAxBvals)\(w*normvdat);
                    end
                    Out.Daxial(xx,yy,zz) = lls(2);
        %            Out.AxSzero(xx,yy,zz) = lls(1);
                    if any(isnan(lls))
                        Out.Daxial(xx,yy,zz) = 0;
       %                 Out.AxSzero(xx,yy,zz) = 0;
                    end
                    
                end

                % radial fitting using single perpendicular direction and
                % monoexponential weighted least squares
                voxdat = squeeze(dwinii.img(xx,yy,zz,RadInd));
                normdat = squeeze(dwinii.img(xx,yy,zz,B0Ind));
                if normdat ~= 0
                    normvdat = log(voxdat./mean(normdat));
                    lls = fitRadBvals\normvdat;
                    % weighted fitting.
                    if ~any(isnan(lls))
                        w = diag(exp(fitRadBvals*lls));
                        lls = (w*fitRadBvals)\(w*normvdat);
                    end
                    Out.Dradial(xx,yy,zz) = lls(2);
         %           Out.RadSzero(xx,yy,zz) = lls(1);
                    if any(isnan(lls))
                        Out.Dradial(xx,yy,zz) = 0;
         %               Out.RadSzero(xx,yy,zz) = 0;
                    end
                end


            end % y loop
        end % x loop
        fprintf('\n');
    end % slice loop

    
    Out.Daxial(Out.Daxial<0) = 0;
    Out.Dradial(Out.Dradial<0) = 0;
    
    % save all of the files in nifti format
    disp('Saving Maps...');

    % This code gets each field from a structure and saves it
    % as a nifti file.  In this way, we don't have to specifically write out each map.  Just 
    % include it as a field and makes saving more versitile under different fitting conditions. 
    niftioutput = fieldnames(Out);
    for i=1:numel(niftioutput)
        niftisave(Out.(niftioutput{i}), OutputPrefix, niftioutput{i}, dwinii);
    end


end %function



function niftisave(imagematrix, Prefix, Suffix, niistruct)
    %function niftisave
    % function that takes in a matrix and saves it to a nii file with
    % proper dimensions.
    %input:
    %  imagematrix - N dimensional matrix (up to 5 below)
    %  filename - Can include path, only nii supported (not nii.gz)
    %  nii - structure from load_nii or similar
    
    fname = strcat(Prefix,'_', Suffix,'.nii');   %Naming convention
    nii = niistruct;
    nii.img=imagematrix;

    % these are essential for double valued output.  scaling factor is removed
    nii.hdr.dime.datatype = 64;
    nii.hdr.dime.bitpix = 64;
    nii.hdr.dime.scl_slope = 1;

    %The following lines edit the hdr for proper sizes
    if size(imagematrix,3)==1
        nii.hdr.dime.dim = [3 size(imagematrix,1) size(imagematrix,2) size(imagematrix,3) 1 1 1 1];
        nii.original.hdr.dime.dim = [3 size(imagematrix,1) size(imagematrix,2) size(imagematrix,3) 1 1 1 1];
    end
    if size(imagematrix,4)==1
        nii.hdr.dime.dim = [3 size(imagematrix,1) size(imagematrix,2) size(imagematrix,3) 1 1 1 1];
        nii.original.hdr.dime.dim = [3 size(imagematrix,1) size(imagematrix,2) size(imagematrix,3) 1 1 1 1];
    end
    if size(imagematrix,4)>1
        nii.hdr.dime.dim = [4 size(imagematrix,1) size(imagematrix,2) size(imagematrix,3) size(imagematrix,4) 1 1 1];
        nii.original.hdr.dime.dim = [4 size(imagematrix,1) size(imagematrix,2) size(imagematrix,3) size(imagematrix,4) 1 1 1];
    end
    if size(imagematrix,5)>1
        nii.hdr.dime.dim = [5 size(imagematrix,1) size(imagematrix,2) size(imagematrix,3) size(imagematrix,4) size(imagematrix,5) 1 1];
        nii.original.hdr.dime.dim = [5 size(imagematrix,1) size(imagematrix,2) size(imagematrix,3)...
            size(imagefile,4) size(imagematrix,5) 1 1];
    end
    fprintf('Writing %s\n',fname);
    save_untouch_nii(nii,fname);


end %function


