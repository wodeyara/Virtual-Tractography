% function that takes as inputs the full interpolated steamlines, any
% parcellation in the MNI 152 space  and an
% optional argument that takes in a lesion mask. In the absence of the
% lesion mask we don't run through the code that checks whether each
% streamline intersects with the lesion mask volume. 
% Output will be the connectome -weighted. Written out in detail below.

function [fullConnectome,redConnectome,noEndStreamlines,eachROIstreamlines, lesMaskStreamlines,ROI2ROI] ...
    = getLesionMaskConnectome(tract_interp, parcels, lesMask)
% INPUTS
% tract_interp: the set of interpolated tracts (say using cubic splines) of
% the form numVertices x 3 (x,y,z) x tracts
% parcels: this is the atlas in the same space as tracts that will be used
% to generate the connectome and the associated outputs, should be a volume
% if read in using read_avw() dont forget to reverse the axes as necessary since it
% comes in reversed for some reason. Compare with tracts and the Brain to
% know which way to rotate and orient the brain.
% lesMask:volume with defined lesion mask where values are all non-zero
% (preferably 1 for ease of understanding) 

% OUTPUTS
% fullConnectome: weighted adjacency matrix of same size as the number of parcels
% redConnectome: weighted adjacency matrix from applying the lesion mask to
% streamlines to identify the SC under injury
% noEndStreamlines: streamlines that don't have termination points inside
% an ROI
% eachROIstreamlines: the set of streamlines emanating from an ROI
% lesMaskStreamlines: the set of streamlines passing through the lesion
% mask
% ROI2ROI: an nx2 matrix where n is number of streamlines, that tells us
% the set of ROIs connected by each streamline

if nargin<2
    disp('Please enter tracts and parcels')
    return
end

if nargin<3
    if length(unique(parcels(:)))>2 && ndims(tract_interp)==3
        lesMask = zeros(size(parcels,1),size(parcels,2) , size(parcels,3));
    else
        disp('Please enter appropriate tracts and parcels')
        return
    end
end

voxel_size = 1; %assuming 1mm x 1mm x 1mm
lengthTract = sum(sqrt(squeeze(sum((tract_interp(2:end,1:3,:) - tract_interp(1:(end-1),1:3,:)).^2, 2))), 1);
disp('assuming tracts are centered at zeros and correcting that so that we can extract volume values')
%if using allTractsInterp this^^ will be true by default
pause(2)
tract_interp(:,1,:) = tract_interp(:,1,:) + 90.5;
tract_interp(:,2,:) = tract_interp(:,2,:) + 108.5;
tract_interp(:,3,:) = tract_interp(:,3,:) + 90.5;

numROIs = length(nonzeros(unique(parcels(:))));
% disp('Number of ROIs is:')
% numROIs
fullConnectome = zeros(numROIs);
redConnectome = zeros(numROIs);

noEndStreamlines = zeros(length(tract_interp),1);

for j =1:numROIs
    eachROIstreamlines{j} = [];
end

lesMaskStreamlines = zeros(size(tract_interp,3),1);

for iTrk=1:size(tract_interp,3)
    
    % Translate continuous vertex coordinates into discrete voxel coordinates
    vox = ceil(squeeze(tract_interp(:,1:3,iTrk)) ./ repmat(voxel_size, size(tract_interp,1),3)); %could multiple tracks wind up with the same voxel coordinate?
    
    % Index into volume to extract scalar values
    inds                = sub2ind([size(parcels,1),size(parcels,2), size(parcels,3)], vox(:,1), vox(:,2), vox(:,3));
    scalars             = parcels(inds);
    % keeping only end points along streamline 
    scalars             = nonzeros([scalars(1),  scalars(100)]);
    if length(unique(scalars))==2
        lesMaskIntersect    = lesMask(inds(1 : 100)); 
    else
        lesMaskIntersect =0;
    end
    % checking for intersections of lesion mask with this streamline anywhere along its path
    
    if length(nonzeros(unique(scalars)))<=1 % checking if some streamlines don't end in two ROIs to form a connection
        noEndStreamlines(iTrk) = 1; 
    else
        % if this works scalars will have the unique ROI numbers that connect
        % any particular streamline. 
        fullConnectome(nonzeros(unique(scalars)),nonzeros(unique(scalars))) =fullConnectome(nonzeros(unique(scalars)),nonzeros(unique(scalars))) + 1/lengthTract(iTrk);% /lengthTract(iTrk) % normalizing by length
        if sum(lesMaskIntersect(:))==0
            redConnectome(nonzeros(unique(scalars)),nonzeros(unique(scalars))) = redConnectome(nonzeros(unique(scalars)),nonzeros(unique(scalars))) + 1/lengthTract(iTrk); %/lengthTract(iTrk)
        else
            lesMaskStreamlines(iTrk) = 1;  % when there is an intersection      
        end
        ROI2ROI(iTrk,1:2) = nonzeros(unique(scalars));
    end
    %% SOME DEBUGGING CODE
    % saving out the streamlines passing through each ROI
    for j = nonzeros(unique(scalars))'
        tmp = eachROIstreamlines{j};
        tmp = [tmp,iTrk];
        eachROIstreamlines{j} = tmp;            
    end
end