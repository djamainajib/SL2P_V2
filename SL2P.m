function SL2P(varargin)

%% 1. Initialization
if ~ismember(nargin,[2,3]), disp({'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!';'--usage : Matlab SL2P_v_01 [input_path\] [S2 tiff data folder] [output_path\ (optional)]'});return; end;

addpath(genpath('.\SL2P_V1-master\tools'));
addpath(genpath('.\SL2P_V1-master\tools\aux_data'));

bio_vars={'LAI','FCOVER','FAPAR','LAI_Cab','LAI_Cw'};
BIO_VAR_bounding_box=importdata('G:\Najib\6_SL2P_regularizing\SL2P_V1-master\tools\aux_data\BIO_VAR_bounding_box.mat');

if nargin==3,   out_path=[varargin{3},strrep(varargin{2}(1:end-5),'L2A','L2B'),'\'];
    else,out_path=[varargin{1},strrep(varargin{2}(1:end-5),'L2A','L2B'),'\'];
end;
if ~isfolder(out_path), mkdir (out_path); end;   
%% 2.1 Loading data........................................................
disp({'===============',varargin{2},'==============='});
disp({'--Loading data--------------------------------------'});
Input_NNT=[]; 
h = waitbar(0,'Loading data...');

L2A_data=[varargin{1},varargin{2}];
L2A_data_name=strsplit(L2A_data,'_');

if sum(ismember(L2A_data_name,{'PRD'}))==1,
   SR=read_S2_MSIL2A_safe_PRD_struct(L2A_data);
else
   SR=read_S2_MSIL2A_safe_struct(L2A_data);
end;
waitbar(1/3)
for bb={'B03','B04','B05','B06','B07','B8A','B11','B12','view_zenith','sun_zenith','view_azimuth','sun_azimuth'}
    eval(['band=SR.',char(bb),';']);
    [r,c]=size(band);
    Input_NNT= [Input_NNT,double(reshape(band,r*c,1))]; 
end;

waitbar(2/3)
%% 2.2 Adding image cordinates
Input_NNT=[reshape((1:r)'*ones(1,c),r*c,1),reshape(ones(1,r)'*(1:c),r*c,1),Input_NNT];
%% 2.3 Organizing input data for NNET (NNET_IN)
Input_NNT(:,end-1)=abs(Input_NNT(:,end-1)-Input_NNT(:,end));Input_NNT(:,end)=[];
Input_NNT(:,3:end-3)=Input_NNT(:,3:end-3)/10000;
Input_NNT(:,end-2:end)=cos(deg2rad(Input_NNT(:,end-2:end))); 
NNT_IN=Input_NNT(:,3:end)';
waitbar(3/3)
close(h)
%% 2.4 Computing input_flags 
input_out_of_range=input_out_of_range_flag_function_SL2P(Input_NNT(:,3:end-3),r,c);% the database is hard coded in the function
%% 2.5 Creating no_bare_soil_or_vegetated_area flag
%% 3. Loading NET
disp({'--Loading NNET--------------------------------------'});
NET_estim=importdata('aux_data\SL2P_NETs.mat');
NET_uncer=importdata('aux_data\SL2P_uncert_NETs.mat');
%% 5. Simulating biophysical parameters (SL2P).....................................
disp({'--Simulating vegetation biophysical variables ------'});
h = waitbar(0,'Simulating bio- variables...');
for ivar=1:length(bio_vars),
    waitbar(ivar/length(bio_vars))
    bio=bio_vars{ivar};
    bio_sim= [Input_NNT(:,1:2),NaN+Input_NNT(:,1)];

    eval(['NET_ivar= NET_estim.',bio,'.NET;']);
    bio_sim (:,3)= sim(NET_ivar, NNT_IN)';
    
    eval(['NET_unc= NET_uncer.',bio,'.NET;']);
    eval(['ps= NET_uncer.',bio,'.Norm_Input;']);
    eval(['ts= NET_uncer.',bio,'.Norm_Output;']);

    NNT_IN_P=mapminmax('apply',NNT_IN,ps);
    bio_sim (:,4)= mapminmax('reverse',sim(NET_unc, NNT_IN_P),ts)';

    %% Creating output_flags
    eval(['bounding_box=BIO_VAR_bounding_box.',bio,';']);
    outpout_out_of_range=outpout_out_of_range_flag_function_SL2P (bio_sim,bounding_box);
    
    bio_sim(find(outpout_out_of_range(:,1)==1),3)=bounding_box.Pmin;
    bio_sim(find(outpout_out_of_range(:,2)==1),3)=bounding_box.Pmax;

    output_thresholded_to_min_outpout= reshape(outpout_out_of_range(:,1),r,c);
    output_thresholded_to_max_outpout= reshape(outpout_out_of_range(:,2),r,c);
    
    output_too_low= reshape(outpout_out_of_range(:,3),r,c);
    output_too_high= reshape(outpout_out_of_range(:,4),r,c);
    %% *********
    flags=(2^0)*input_out_of_range+(2^1)*output_thresholded_to_min_outpout+(2^2)*output_thresholded_to_max_outpout+(2^3)*output_too_low+(2^4)*output_too_high;

    NNT_OUT=[];
    NNT_OUT.xb=SR.xb;
    NNT_OUT.yb=SR.yb;
    NNT_OUT.lat=SR.lat;
    NNT_OUT.lon=SR.lon;
    NNT_OUT.Ib=SR.Ib;

    eval(['NNT_OUT.',lower(bio),'=reshape(bio_sim(:,3),r,c);']);
    eval(['NNT_OUT.',lower(bio),'_Uncertainties=reshape(bio_sim(:,4),r,c);']);
    eval(['NNT_OUT.',lower(bio),'_flags=flags;']);
    
    eval(['NNT_OUT.',lower(bio),'_input_out_of_range= input_out_of_range;']);
    eval(['NNT_OUT.',lower(bio),'_output_thresholded_to_min_outpout= output_thresholded_to_min_outpout;']);
    eval(['NNT_OUT.',lower(bio),'_output_thresholded_to_max_outpout= output_thresholded_to_max_outpout;']);
    eval(['NNT_OUT.',lower(bio),'_output_too_low= output_too_low;']);
    eval(['NNT_OUT.',lower(bio),'_output_too_high= output_too_high;']);
    %% exporting tif files
    %bbox=Ib.BoundingBox;
    %[bbox(:,2),bbox(:,1)] = utm2deg(bbox(:,1),bbox(:,2),repmat(NNT_OUT.Ib.utmzone,2,1));    
    bit_depth=32;
    geotiffwrite([out_path,strrep(varargin{2}(1:end-5),'L2A','L2B'),'_',lower(bio),'.tif'], bbox, eval(['NNT_OUT.',lower(bio)]), bit_depth, Ib);
    geotiffwrite([out_path,strrep(varargin{2}(1:end-5),'L2A','L2B'),'_',lower(bio),'_uncertainties.tif'], bbox, eval(['NNT_OUT.',lower(bio),'_Uncertainties']), bit_depth, Ib);
    geotiffwrite([out_path,strrep(varargin{2}(1:end-5),'L2A','L2B'),'_',lower(bio),'_flags.tif'], bbox, eval(['NNT_OUT.',lower(bio),'_flags']), bit_depth, Ib);  
    % exporting .mat files
    save([out_path,strrep(varargin{2}(1:end-5),'L2A','L2B'),'_',lower(bio),'.mat'],'NNT_OUT','-v7.3');
end;
geotiffwrite([out_path,strrep(varargin{2}(1:end-5),'L2A','L2B'),'_lat_map.tif'], bbox, eval(['NNT_OUT.lat']), bit_depth, NNT_OUT.Ib);
geotiffwrite([out_path,strrep(varargin{2}(1:end-5),'L2A','L2B'),'_lon_map.tif'], bbox, eval(['NNT_OUT.lon']), bit_depth, NNT_OUT.Ib);  

geotiffwrite([out_path,'lat_map.tif'], bbox, eval(['NNT_OUT.lat']), bit_depth, NNT_OUT.Ib);
geotiffwrite([out_path,'lon_map.tif'], bbox, eval(['NNT_OUT.lon']), bit_depth, NNT_OUT.Ib);  
    
disp({'--End SL2P ------'});
close(h)
end


