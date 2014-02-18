%%readhe5_sp_nwus_wcld
%%arr 04/28/2011

%run in command prompt using batch file read_nwus.bat 
%(C:\Program Files\MATLAB\R2008a) to avoid memory issues

addpath C:\Ashley\tools\Regrid_tools
addpath C:\Ashley\tools\Regrid_tools\Compute_Corner_Pts_ms


warning off all
tic

%****************************%
lonmin = -125;  lonmax = -95;
latmin = 37.5;    latmax = 50;
%****************************%

satellite='OMI';
retrieval='SP';

Latitude=zeros(60,2000);
Longitude=zeros(60,300);
SpacecraftAltitude=zeros(300,1);
SpacecraftLatitude=zeros(300,1);
SpacecraftLongitude=zeros(300,1);
Time=zeros(300,1);
ViewingZenithAngle=zeros(60,300);
SolarZenithAngle=zeros(60,300);
ViewingAzimuthAngle=zeros(60,300);
SolarAzimuthAngle=zeros(60,300);
AMFInitial=zeros(60,300);
AMFPolluted=zeros(60,300);
AMFUnpolluted=zeros(60,300);
CloudFraction=zeros(60,300);
CloudRadianceFraction=zeros(60,300);
ColumnAmountNO2=zeros(60,300);
ColumnAmountNO2Initial=zeros(60,300);
ColumnAmountNO2Polluted=zeros(60,300);
SlantColumnAmountNO2=zeros(60,300);
TerrainHeight=zeros(60,300);
TerrainPressure=zeros(60,300);
TerrainReflectivity=zeros(60,300);
vcdQualityFlags=zeros(60,300);
CloudPressure=zeros(60,300);
ColumnAmountNO2Trop=zeros(60,300);
TropFractionUnpolluted=zeros(60,300);

date_start='2004/10/01';
date_end='2011/09/29';

last_file=dir(fullfile('J:','New_Retrieval_sp_mats','NWUS_wcld'));
if length(last_file(end-1).name)<5;
    last_date='2004/09/30';
else
    last_date=last_file(end).name(8:15); last_date=([last_date(1:4),'/',last_date(5:6),'/',last_date(7:8)]);
end

total_days=datenum(date_end)-datenum(last_date)+1;
%last_date='2007/12/18';
for j=1:total_days;
    R=addtodate(datenum(last_date), j, 'day');
    date=datestr(R,26);
    year=date(1:4);
    month=date(6:7);
    day=date(9:10);
    
    directory=(['I:\OMI_SP_NO2\',year,'\',month]);
    cd(directory)    
    Data=struct('Longitude',0,'Latitude',0,'Loncorn',0,'Latcorn',0,'Time',0,'ViewingZenithAngle',0,'SolarZenithAngle',0,'ViewingAzimuthAngle',0,'SolarAzimuthAngle',0,'AMFPolluted',0,'AMFUnpolluted',0,'CloudFraction',0,'CloudRadianceFraction',0,'ColumnAmountNO2Polluted',0,'TerrainHeight',0,'TerrainPressure',0,'TerrainReflectivity',0,'vcdQualityFlags',0,'CloudPressure',0);
    file=['OMI-Aura_L2-OMNO2_',year,'m',month,day,'*.he5']; 
    sp_files=dir(fullfile('I:','OMI_SP_NO2',year,month,file));
    nn=length(sp_files);
    sp_files(nn+1).name='OMI-Aura_L2-OMNO2_0000m0000t9999-o00000_v000-0000m0000t000000.he5';
    a=sp_files;
    cd(directory)  
    E=0;
    if isempty(a);
        disp(['No Data Available For ',month,' ',day,' ',year])
    else
        for e=1:nn    %for loop over all swaths from a day 
            filename=sp_files(e).name;
            %if str2double(filename(29:30))<18||str2double(filename(29:30))>22%||exist(filename2,'file')==0;
                %continue
            %else      
            cd(directory)
            hinfo=hdf5info(filename);
                
            %Latitude
            Latitude=hdf5read(hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(2).Datasets(2)); 
            Row=0:59; Row=Row'; Row=repmat(Row,1,size(Latitude,2));
            Swath=filename(35:39); Swath=str2double(Swath).*ones(size(Latitude));
            
            lat=Latitude';
            lat_i=[37.5 50];
            [i_i j_j]=find(lat > lat_i(1) - 0.25 & lat < lat_i(2) + 0.25);
            cut_y=min(i_i):max(i_i);
            cut_x = 1:60;
            lat=double(lat(cut_y,cut_x));
            Latitude=Latitude(cut_x,cut_y)'; Latitude=double(Latitude); 
            Row=Row(cut_x,cut_y)';
            Swath=Swath(cut_x,cut_y)';
                  
            stride = [1 1 1];
            blocksize = [1 1 1];
            offset = [(min(i_i)-1),0]; 
            slabsize = [length(cut_y),60];
            memspaceID = H5S.create_simple(length(slabsize), slabsize, slabsize);
            fileID = H5F.open(filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
                   
            %Longitude
            datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(2).Datasets(3).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); Longitude = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); Longitude=double(Longitude); Longitude=Longitude'; 
            %ViewingAzimuthAngle
            datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(2).Datasets(10).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); ViewingAzimuthAngle = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); ViewingAzimuthAngle=double(ViewingAzimuthAngle); ViewingAzimuthAngle=ViewingAzimuthAngle';
            %ViewingZenithAngle
            datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(2).Datasets(11).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); ViewingZenithAngle = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); ViewingZenithAngle=double(ViewingZenithAngle); ViewingZenithAngle=ViewingZenithAngle'; 
            %SolarAzimuthAngle
            datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(2).Datasets(4).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); SolarAzimuthAngle = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); SolarAzimuthAngle=double(SolarAzimuthAngle); SolarAzimuthAngle=SolarAzimuthAngle';                 
            %SolarZenithAngle
            datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(2).Datasets(5).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); SolarZenithAngle = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); SolarZenithAngle=double(SolarZenithAngle); SolarZenithAngle=SolarZenithAngle';  
                   
            slabsize = [length(cut_y),1];
            memspaceID = H5S.create_simple(length(slabsize), slabsize, slabsize);
                    
                    
            %SpacecraftAltitude
            datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(2).Datasets(6).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); SpacecraftAltitude = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); SpacecraftAltitude=double(SpacecraftAltitude); SpacecraftAltitude=SpacecraftAltitude';   
            %SpacecraftLatitude
            datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(2).Datasets(7).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); SpacecraftLatitude = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); SpacecraftLatitude=double(SpacecraftLatitude); SpacecraftLatitude=SpacecraftLatitude';
            %SpacecraftLongitude
            datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(2).Datasets(8).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); SpacecraftLongitude = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); SpacecraftLongitude=double(SpacecraftLongitude); SpacecraftLongitude=SpacecraftLongitude';
            %Time
            datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(2).Datasets(9).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); Time = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); Time=double(Time); Time=Time'; 
            Time=repmat(Time,1,60);

            lat=Latitude; 
            lon=Longitude;
            Lat=lat; Lon=Longitude;
            x=find(Lon>lonmax | Lon<lonmin);
            y=find(Lat>latmax | Lat<latmin);
            Lon(x)=NaN;     Lon(y)=NaN;     Lon(isnan(Lon))=[]; 
            Lat(x)=NaN;     Lat(y)=NaN;     Lat(isnan(Lat))=[];
    
            if isempty(Lon)==1 || isempty(Lat)==1 || length(Lat)==1; 
                continue
            else
                corner_coordinates
                E=E+1;
                lat = corners(:,:,2,5);
                lon = corners(:,:,1,5);
                latcorn = corners(:,:,2,1:4); latcorn = squeeze(latcorn); 
                a = latcorn(:,:,1); a = a(:); b = latcorn(:,:,2); b = b(:); c = latcorn(:,:,3); c = c(:); d = latcorn(:,:,4); d = d(:);
                latcorn = [a,b,c,d]; latcorn = latcorn';
                loncorn = corners(:,:,1,1:4); loncorn = squeeze(loncorn);
                a = loncorn(:,:,1); a = a(:); b = loncorn(:,:,2); b = b(:); c = loncorn(:,:,3); c = c(:); d = loncorn(:,:,4); d = d(:);
                loncorn = [a,b,c,d]; loncorn = loncorn';
                        
                slabsize = [length(cut_y),60];
                memspaceID = H5S.create_simple(length(slabsize), slabsize, slabsize);
                      
                %AMFInitial
                datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(1).Datasets(1).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); AMFInitial = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); AMFInitial=double(AMFInitial); AMFInitial=AMFInitial';
                %AMFPolluted
                datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(1).Datasets(7).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); AMFPolluted = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); AMFPolluted=double(AMFPolluted); AMFPolluted=AMFPolluted';
                %AMFUnpolluted
                datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(1).Datasets(16).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); AMFUnpolluted = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); AMFUnpolluted=double(AMFUnpolluted); AMFUnpolluted=AMFUnpolluted';
                %CloudFraction
                datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(1).Datasets(23).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); CloudFraction = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); CloudFraction=double(CloudFraction); CloudFraction=CloudFraction';
                %CloudPressure
                datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(1).Datasets(25).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); CloudPressure = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); CloudPressure=double(CloudPressure); CloudPressure=CloudPressure';
                %CloudRadianceFraction
                datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(1).Datasets(27).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); CloudRadianceFraction = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); CloudRadianceFraction=double(CloudRadianceFraction); CloudRadianceFraction=CloudRadianceFraction';
                %ColumnAmountNO2
                datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(1).Datasets(28).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); ColumnAmountNO2 = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); ColumnAmountNO2=double(ColumnAmountNO2); ColumnAmountNO2=ColumnAmountNO2';
                %ColumnAmountNO2Initial
                datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(1).Datasets(31).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); ColumnAmountNO2Initial = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); ColumnAmountNO2Initial=double(ColumnAmountNO2Initial); ColumnAmountNO2Initial=ColumnAmountNO2Initial';
                %ColumnAmountNO2Polluted
                datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(1).Datasets(33).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); ColumnAmountNO2Polluted = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); ColumnAmountNO2Polluted=double(ColumnAmountNO2Polluted); ColumnAmountNO2Polluted=ColumnAmountNO2Polluted';
                %SlantColumnAmountNO2
                datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(1).Datasets(50).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); SlantColumnAmountNO2 = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); SlantColumnAmountNO2=double(SlantColumnAmountNO2); SlantColumnAmountNO2=SlantColumnAmountNO2';
                %TerrainHeight
                datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(1).Datasets(58).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); TerrainHeight = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); TerrainHeight=double(TerrainHeight); TerrainHeight=TerrainHeight';
                %TerrainPressure
                datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(1).Datasets(59).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); TerrainPressure = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); TerrainPressure=double(TerrainPressure); TerrainPressure=TerrainPressure';
                %TerrainReflectivity
                datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(1).Datasets(60).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); TerrainReflectivity = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); TerrainReflectivity=double(TerrainReflectivity); TerrainReflectivity=TerrainReflectivity';
               if str2double(filename(35:39))<24349;
                    %vcdQualityFlags
                    datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(1).Datasets(67).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); vcdQualityFlags = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); vcdQualityFlags=double(vcdQualityFlags); vcdQualityFlags=vcdQualityFlags';
                else
                    %vcdQualityFlags
                    datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(1).Datasets(68).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); vcdQualityFlags = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); vcdQualityFlags=double(vcdQualityFlags); vcdQualityFlags=vcdQualityFlags';
               end
               %ColumnAmountNO2Trop
                datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(1).Datasets(36).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); ColumnAmountNO2Trop = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); ColumnAmountNO2Trop=double(ColumnAmountNO2Trop); ColumnAmountNO2Trop=ColumnAmountNO2Trop';
                %TropFractionUnpolluted
                datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(1).Datasets(61).Name); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); TropFractionUnpolluted = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); TropFractionUnpolluted=double(TropFractionUnpolluted); TropFractionUnpolluted=TropFractionUnpolluted';

                H5F.close(fileID); %close omi file to free up space

                s=size(SolarAzimuthAngle);
                latcorn=reshape(latcorn,4,s(1),s(2));
                loncorn=reshape(loncorn,4,s(1),s(2));
                Time=reshape(Time,s(1),s(2));
                        

                %RelativeAzimuthAngle=abs(SolarAzimuthAngle-ViewingAzimuthAngle);
                %i_raa=find(RelativeAzimuthAngle > 180);
                %RelativeAzimuthAngle(i_raa)=360-RelativeAzimuthAngle(i_raa);
                RelativeAzimuthAngle=abs(SolarAzimuthAngle+180-ViewingAzimuthAngle);
                RelativeAzimuthAngle(RelativeAzimuthAngle > 180)=360-RelativeAzimuthAngle(RelativeAzimuthAngle > 180);
                
                        
                pol_i=find(ColumnAmountNO2Polluted<-1E16);
                ColumnAmountNO2Polluted(pol_i)=0;
                loncorn=loncorn(1:4,:);
                latcorn=latcorn(1:4,:);

                         
                x=find(lon<lonmin | lon>lonmax);
                y=find(lat<latmin | lat>latmax);
                lon(x)=NaN;                         lon(y)=NaN;                         lon(isnan(lon))=[]; 
                lat(x)=NaN;                         lat(y)=NaN;                         lat(isnan(lat))=[];
                loncorn(:,x)=NaN;                   loncorn(:,y)=NaN;                   loncorn=loncorn';                   loncorn(any(isnan(loncorn)'),:) = [];               loncorn=loncorn';
                latcorn(:,x)=NaN;                   latcorn(:,y)=NaN;                   latcorn=latcorn';                   latcorn(any(isnan(latcorn)'),:) = [];               latcorn=latcorn';
                SolarAzimuthAngle(x)=NaN;           SolarAzimuthAngle(y)=NaN;           SolarAzimuthAngle(isnan(SolarAzimuthAngle))=[];
                SolarZenithAngle(x)=NaN;            SolarZenithAngle(y)=NaN;            SolarZenithAngle(isnan(SolarZenithAngle))=[];
                ViewingAzimuthAngle(x)=NaN;         ViewingAzimuthAngle(y)=NaN;         ViewingAzimuthAngle(isnan(ViewingAzimuthAngle))=[];
                ViewingZenithAngle(x)=NaN;          ViewingZenithAngle(y)=NaN;          ViewingZenithAngle(isnan(ViewingZenithAngle))=[];
                Time(x)=NaN;                        Time(y)=NaN;                        Time(isnan(Time))=[];
                AMFInitial(x)=NaN;                  AMFInitial(y)=NaN;                  AMFInitial(isnan(AMFInitial))=[];
                AMFPolluted(x)=NaN;                 AMFPolluted(y)=NaN;                 AMFPolluted(isnan(AMFPolluted))=[];
                AMFUnpolluted(x)=NaN;               AMFUnpolluted(y)=NaN;               AMFUnpolluted(isnan(AMFUnpolluted))=[];
                CloudFraction(x)=NaN;               CloudFraction(y)=NaN;               CloudFraction(isnan(CloudFraction))=[];
                CloudPressure(x)=NaN;               CloudPressure(y)=NaN;               CloudPressure(isnan(CloudPressure))=[];
                CloudRadianceFraction(x)=NaN;       CloudRadianceFraction(y)=NaN;       CloudRadianceFraction(isnan(CloudRadianceFraction))=[];
                ColumnAmountNO2(x)=NaN;             ColumnAmountNO2(y)=NaN;             ColumnAmountNO2(isnan(ColumnAmountNO2))=[];
                ColumnAmountNO2Initial(x)=NaN;      ColumnAmountNO2Initial(y)=NaN;      ColumnAmountNO2Initial(isnan(ColumnAmountNO2Initial))=[];
                ColumnAmountNO2Polluted(x)=NaN;     ColumnAmountNO2Polluted(y)=NaN;     ColumnAmountNO2Polluted(isnan(ColumnAmountNO2Polluted))=[];
                SlantColumnAmountNO2(x)=NaN;        SlantColumnAmountNO2(y)=NaN;        SlantColumnAmountNO2(isnan(SlantColumnAmountNO2))=[];
                TerrainHeight(x)=NaN;               TerrainHeight(y)=NaN;               TerrainHeight(isnan(TerrainHeight))=[];
                TerrainPressure(x)=NaN;             TerrainPressure(y)=NaN;             TerrainPressure(isnan(TerrainPressure))=[];
                TerrainReflectivity(x)=NaN;         TerrainReflectivity(y)=NaN;         TerrainReflectivity(isnan(TerrainReflectivity))=[];
                vcdQualityFlags(x)=NaN;             vcdQualityFlags(y)=NaN;             vcdQualityFlags(isnan(vcdQualityFlags))=[];
                RelativeAzimuthAngle(x)=NaN;        RelativeAzimuthAngle(y)=NaN;        RelativeAzimuthAngle(isnan(RelativeAzimuthAngle))=[];
                ColumnAmountNO2Trop(x)=NaN;         ColumnAmountNO2Trop(y)=NaN;         ColumnAmountNO2Trop(isnan(ColumnAmountNO2Trop))=[];
                TropFractionUnpolluted(x)=NaN;      TropFractionUnpolluted(y)=NaN;      TropFractionUnpolluted(isnan(TropFractionUnpolluted))=[];
                Row(x)=NaN;                         Row(y)=NaN;                         Row(isnan(Row))=[];                      
                Swath(x)=NaN;                       Swath(y)=NaN;                       Swath(isnan(Swath))=[];
                           

                Data(E).Latitude = lat(:);                                  Data(E).AMFPolluted = AMFPolluted(:);
                Data(E).Longitude = lon(:);                                 Data(E).AMFUnpolluted = AMFUnpolluted(:);
                Data(E).Loncorn = loncorn(1:4,:);                             Data(E).CloudFraction = CloudFraction(:);
                Data(E).Latcorn = latcorn(1:4,:);                             Data(E).ColumnAmountNO2Polluted = ColumnAmountNO2Polluted(:);
                Data(E).SolarAzimuthAngle = SolarAzimuthAngle(:);                   
                Data(E).SolarZenithAngle = SolarZenithAngle(:);             Data(E).SlantColumnAmountNO2 = SlantColumnAmountNO2(:);
                Data(E).ViewingAzimuthAngle = ViewingAzimuthAngle(:);       Data(E).TerrainHeight = TerrainHeight(:);
                Data(E).ViewingZenithAngle = ViewingZenithAngle(:);         Data(E).TerrainPressure = TerrainPressure(:);
                Data(E).Time = Time(:);                                     Data(E).TerrainReflectivity = TerrainReflectivity(:);
                Data(E).ColumnAmountNO2 = ColumnAmountNO2(:);               Data(E).vcdQualityFlags = vcdQualityFlags(:);
                Data(E).AMFInitial = AMFInitial(:);                         Data(E).CloudRadianceFraction = CloudRadianceFraction(:);
                Data(E).ColumnAmountNO2Initial = ColumnAmountNO2Initial(:); Data(E).RelativeAzimuthAngle = RelativeAzimuthAngle(:);
                Data(E).CloudPressure = CloudPressure(:);                   Data(E).ColumnAmountNO2Trop = ColumnAmountNO2Trop(:);
                Data(E).TropFractionUnpolluted = TropFractionUnpolluted(:);
                Data(E).Row = Row(:);
                Data(E).Swath = Swath(:);
                        
 %Add MODIS cloud info to the files%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       
                        
                d2=1+datenum(str2double(year),str2double(month),str2double(day))-datenum(str2double(year),1,1);
                x=numel(num2str(d2));
                if x==1;
                    day2=(['00',num2str(d2)]);
                elseif x==2;
                    day2=(['0',num2str(d2)]);
                elseif x==3;
                    day2=num2str(d2);
                end    
                modis_file=(['MYD06_L2.A',year,day2,'*.hdf']);
                %modis_files=dir(fullfile('\\128.32.208.19\Elements (G)\MODIS_Cloud_US',year,modis_file));
                modis_files=dir(fullfile('I:\MODIS_Cloud_US',year,modis_file));
                n=length(modis_files);
                %mod_directory=['\\128.32.208.19\Elements (G)\MODIS_Cloud_US\',year];
                mod_directory=['I:\MODIS_Cloud_US\',year];
                cd(mod_directory) 
                for ii=1:n;
                    mod_filename=modis_files(ii).name;
                    if str2double(mod_filename(19:22))<str2double(sp_files(e).name(29:32));
                       continue
                    elseif str2double(mod_filename(19:22))>str2double(sp_files(e+1).name(29:32));
                       continue
                    else
                        mod_fileinfo=hdfinfo(mod_filename);
                        Latitude=hdfread(mod_fileinfo.Vgroup(1).Vgroup(1).SDS(1)); Latitude=double(Latitude); Latitude=Latitude(:);
                        Longitude=hdfread(mod_fileinfo.Vgroup(1).Vgroup(1).SDS(2)); Longitude=double(Longitude); Longitude=Longitude(:);
                        CloudFraction=hdfread(mod_fileinfo.Vgroup(1).Vgroup(2).SDS(18)); CloudFraction=double(CloudFraction); CloudFraction=CloudFraction(:); 
                        CloudFraction(CloudFraction==127)=100; CloudFraction=CloudFraction*0.009999999776482582;

                        x=find(Longitude>lonmax | Longitude<lonmin);
                        y=find(Latitude>latmax | Latitude<latmin);
                        Longitude(x)=NaN;           Longitude(y)=NaN;               Longitude(isnan(Longitude))=[]; 
                        Latitude(x)=NaN;            Latitude(y)=NaN;                Latitude(isnan(Latitude))=[];
                        CloudFraction(x)=NaN;       CloudFraction(y)=NaN;           CloudFraction(isnan(CloudFraction))=[];

                        if isempty(Longitude)||isempty(Latitude);
                        else
                            if exist('mod_Data','var')==0;
                                mod_Data(1).Longitude=Longitude;              mod_Data(1).Latitude=Latitude;    
                                mod_Data(1).CloudFraction=CloudFraction;     
                            elseif exist('mod_Data','var')==1;
                                mod_Data(1).Longitude=[mod_Data(1).Longitude;Longitude];
                                mod_Data(1).Latitude=[mod_Data(1).Latitude;Latitude];
                                mod_Data(1).CloudFraction=[mod_Data(1).CloudFraction;CloudFraction];
                            end
                        end
                    end
                end
                %if exist('mod_Data','var')==0
                %else
                %    X=24.025:0.05:50.975;
                %    Y=-125.975:0.05:-64.025;
                %    [YI,XI]=meshgrid(Y,X);
                %    ZI = griddata(mod_Data.Longitude,mod_Data.Latitude,mod_Data.CloudFraction,YI,XI);
                %    mod_Data.CloudFractionGrid=ZI;
                %end
                       
                if exist('mod_Data','var')==0;
                    Data(E).NewCloud=-127*ones(length(Data(E).Latitude),1);
                else
                    for jj=1:length(Data(E).Latitude);
                        x = [];
                        x1 = Data(E).Loncorn(1,jj);   y1 = Data(E).Latcorn(1,jj);
                        x2 = Data(E).Loncorn(2,jj);   y2 = Data(E).Latcorn(2,jj);
                        x3 = Data(E).Loncorn(3,jj);   y3 = Data(E).Latcorn(3,jj);
                        x4 = Data(E).Loncorn(4,jj);   y4 = Data(E).Latcorn(4,jj);

                        xall=[x1;x2;x3;x4;x1];
                        yall=[y1;y2;y3;y4;y1];
                        %get modis Cld
                        xx = inpolygon(mod_Data.Latitude,mod_Data.Longitude,yall,xall);

                        cld_vals=mod_Data.CloudFraction(xx); %cld_zeros=find(cld_vals==0); %cld_vals(cld_zeros)=NaN; 
                        cld_vals(isnan(cld_vals))=[];
                        Data(E).NewCloud(jj,1)=mean(cld_vals); %Data(E).NewCloud=Data(E).NewCloud';

                        clear lat1 lat2 lat3 lat4 xx
                    end
                end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                    
            end
        end
        cd J:\New_Retrieval_sp_mats\NWUS_wCld
        savename=[satellite,'_',retrieval,'_',year,month,day];
        save(savename, 'Data')
        clear Data 
        pack
        toc
        t=toc;
        if t>1200
            quit
        end
    end
end

