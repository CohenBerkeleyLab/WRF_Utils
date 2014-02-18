%%readhe5_sp_ca
%%arr 11/20/2009

addpath C:\Ashley\NewCaRetrieval2\Regrid_tools
addpath C:\Ashley\NewCaRetrieval2\Regrid_tools\Compute_Corner_Pts_ms


warning off all
tic

%****************************%
lonmin = -126;  lonmax = -113;
latmin = 31;    latmax = 44;
%****************************%

satellite='OMI';
retrieval='SP';

%years={'2004';'2004';'2004';'2005';'2005';'2005';'2005';'2005';'2005';'2005';'2005';'2005';'2005';'2005';'2005';'2006';'2006';'2006';'2006';'2006';'2006';'2006';'2006';'2006';'2006';'2006';'2006';'2007';'2007';'2007';'2007';'2007';'2007';'2007';'2007';'2007';'2007';'2007';'2007';'2008';'2008';'2008';'2008';'2008';'2008';'2008';'2008';'2008';'2008';'2008';'2008';'2009';'2009';'2009';'2009';'2009';'2009';'2009';'2009';'2009';'2009';'2009';'2009';'2010';'2010';'2010';'2010';'2010';'2010'};
%months={'10';'11';'12';'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'01';'02';'03';'04';'05';'06'};
%years={'2005';'2005';'2005';'2005';'2005';'2005';'2005';'2005';'2005';'2005';'2005';'2005';'2006';'2006';'2006';'2006';'2006';'2006';'2006';'2006';'2006';'2006';'2006';'2006';'2007';'2007';'2007';'2007';'2007';'2007';'2007';'2007';'2007';'2007';'2007';'2007';'2008';'2008';'2008';'2008';'2008';'2008';'2008';'2008';'2008';'2008';'2008';'2008'};
%months={'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12'};
%years={'2004';'2004';'2004';'2009';'2009';'2009';'2009';'2009';'2009';'2009';'2009';'2009';'2009';'2009';'2009';'2010';'2010';'2010';'2010';'2010';'2010';'2010';'2010';'2010';'2010';'2010'};
%months={'10';'11';'12';'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11'};
years={'2010'};
months={'03'};

%{
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
%}
for j=1:length(months);
    year=years{j};
    month=months{j};
    eom=eomday(str2double(year),str2double(month));
    days={'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31'};
    directory=(['I:\OMI_SP_NO2\',year,'\',month]);
    cd(directory)    
    for i=6;%:eom;
        day=days{i};
        Data=struct('Longitude',0,'Latitude',0,'Loncorn',0,'Latcorn',0,'Time',0,'ViewingZenithAngle',0,'SolarZenithAngle',0,'ViewingAzimuthAngle',0,'SolarAzimuthAngle',0,'AMFPolluted',0,'AMFUnpolluted',0,'CloudFraction',0,'CloudRadianceFraction',0,'ColumnAmountNO2Polluted',0,'TerrainHeight',0,'TerrainPressure',0,'TerrainReflectivity',0,'vcdQualityFlags',0,'CloudPressure',0);
        file=['OMI-Aura_L2-OMNO2_',year,'m',month,day,'*.he5']; 
        sp_files=dir(fullfile('I:','OMI_SP_NO2',year,month,file));
        n=length(sp_files);
        a=sp_files;
        cd(directory)  
        E=0;
        if isempty(a);
            disp(['No Data Available For ',month,' ',day,' ',year])
        else
            for e=1:n    %for loop over all swaths from a day 
                filename=sp_files(e).name;
                %filename2=['doas_cld_',filename(19:32),'.mat'];
                %directory2=['J:\OMI_MATLAB_CLDS_DOAS\doas_cld\',year,'\',month,'\'];
                %cd(directory2) 
                if str2double(filename(29:30))<18||str2double(filename(29:30))>22%||exist(filename2,'file')==0;
                    continue
                else      
                    cd(directory)
                    hinfo=hdf5info(filename);
                    
                     %Latitude
                    Latitude=hdf5read(hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(2).Datasets(2)); %JLL If this is supposed to be the center latitude, then it needs to be Datasets(4)
                    Row=0:59; Row=Row'; Row=repmat(Row,1,size(Latitude,2));
                    Swath=filename(35:39); Swath=str2double(Swath).*ones(size(Latitude));
            
                    lat=Latitude';
                    lat_i=[30.5 44.5];
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
                    
                    %{
                    l=length(Latitude);
                    A=zeros(l,1);
                    for l_index=1:l;
                        L=find(Latitude(1:60,l_index)>latmin&Latitude(1:60,l_index)<latmax);
                        A(l_index)=sum(L);
                    end
                    B=find(A~=0); C=length(B);
                    Latitude=Latitude(1:60,B(1):B(C)); Latitude=Latitude';
              
                    stride = [1 1 1];
                    blocksize = [1 1 1];
                    offset = [(B(1)-1),0]; 
                    slabsize = [length(B(1):B(C)),60];
                    memspaceID = H5S.create_simple(length(slabsize), slabsize, slabsize);
                    fileID = H5F.open(filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
                    %}
                    
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
                        
                        s=size(SolarAzimuthAngle);
                        latcorn=reshape(latcorn,4,s(1),s(2));
                        loncorn=reshape(loncorn,4,s(1),s(2));
                        Time=reshape(Time,s(1),s(2));
                        
                        %filename2=['doas_cld_',filename(19:32)];
                        %directory2=['J:\OMI_MATLAB_CLDS_DOAS\doas_cld\',year,'\',month,'\'];
                        %cd(directory2) 
                        %%if exist(filename2,'file')==0;
                            %%continue
                        %%else
                        %load(filename2,'N_s_OMI','R_c_OMI','rel_vza')
                        %SlantCol=N_s_OMI(1:size(rel_vza,1),1:size(rel_vza,2));
                        %Radiance=R_c_OMI(1:size(rel_vza,1),1:size(rel_vza,2));
                        SlantCol=zeros(size(ColumnAmountNO2Trop));
                        Radiance=zeros(size(ColumnAmountNO2Trop));

%                         RelativeAzimuthAngle=abs(SolarAzimuthAngle-ViewingAzimuthAngle);
%                         i_raa=find(RelativeAzimuthAngle > 180);
%                         RelativeAzimuthAngle(i_raa)=360-RelativeAzimuthAngle(i_raa);
                        RelativeAzimuthAngle=abs(SolarAzimuthAngle+180-ViewingAzimuthAngle);
                        RelativeAzimuthAngle(RelativeAzimuthAngle > 180)=360-RelativeAzimuthAngle(RelativeAzimuthAngle > 180);
                

                        %{
                        lat_i=[30.5 44.5];
                        [ii jj]=find(lat > lat_i(1) - 0.257 & lat < lat_i(2) + 0.257);
                        cut_y=min(ii):max(ii);
                        cut_x = 1:60;
                        lat=double(lat(cut_y,cut_x));
                        
                        if numel(lat)~=numel(Radiance);
                            lat=Latitude;
                            lat_i=[30.5 44.5];
                            [ii jj]=find(lat > lat_i(1) - 0.25 & lat < lat_i(2) + 0.25);
                            cut_y=min(ii):max(ii);
                            cut_x = 1:60;
                            lat=double(lat(cut_y,cut_x));
                        end
                        %}      
                        %{
                        lon=double(lon(cut_y,cut_x));
                        loncorn=double(loncorn(:,cut_y,cut_x));
                        latcorn=double(latcorn(:,cut_y,cut_x));
                        SolarAzimuthAngle=double(SolarAzimuthAngle(cut_y,cut_x));
                        SolarZenithAngle=double(SolarZenithAngle(cut_y,cut_x));
                        ViewingAzimuthAngle=double(ViewingAzimuthAngle(cut_y,cut_x));
                        ViewingZenithAngle=double(ViewingZenithAngle(cut_y,cut_x));
                        Time=double(Time(cut_y,cut_x));
                        ColumnAmountNO2=double(ColumnAmountNO2(cut_y,cut_x));
                        AMFInitial=double(AMFInitial(cut_y,cut_x));
                        ColumnAmountNO2Initial=double(ColumnAmountNO2Initial(cut_y,cut_x));
                        CloudPressure=double(CloudPressure(cut_y,cut_x));
                        TropFractionUnpolluted=double(TropFractionUnpolluted(cut_y,cut_x));
                        AMFPolluted=double(AMFPolluted(cut_y,cut_x));
                        AMFUnpolluted=double(AMFUnpolluted(cut_y,cut_x));
                        CloudFraction=double(CloudFraction(cut_y,cut_x));
                        ColumnAmountNO2Polluted=double(ColumnAmountNO2Polluted(cut_y,cut_x));
                        SlantColumnAmountNO2=double(SlantColumnAmountNO2(cut_y,cut_x));
                        TerrainHeight=double(TerrainHeight(cut_y,cut_x));
                        TerrainPressure=double(TerrainPressure(cut_y,cut_x));
                        TerrainReflectivity=double(TerrainReflectivity(cut_y,cut_x));
                        vcdQualityFlags=double(vcdQualityFlags(cut_y,cut_x));
                        CloudRadianceFraction=double(CloudRadianceFraction(cut_y,cut_x));
                        RelativeAzimuthAngle=double(RelativeAzimuthAngle(cut_y,cut_x));
                        ColumnAmountNO2Trop=double(ColumnAmountNO2Trop(cut_y,cut_x));
                        %}
                        
                        pol_i=find(ColumnAmountNO2Polluted<-1E16);
                        ColumnAmountNO2Polluted(pol_i)=0;
                        loncorn=loncorn(1:4,:);
                        latcorn=latcorn(1:4,:);
                        %lat(lat<-1E30)=NaN; a=find(lat==0); lat(a)=NaN;
                        %lon(lon<-1E30)=NaN; lon(a)=NaN;
                         
                        x=find(lon<-126 | lon>-113);
                        y=find(lat<31 | lat>44);
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
                        SlantCol(x)=NaN;                    SlantCol(y)=NaN;                    SlantCol(isnan(SlantCol))=[];
                        Radiance(x)=NaN;                    Radiance(y)=NaN;                    Radiance(isnan(Radiance))=[];
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
                        Data(E).SlantCol = SlantCol(:);
                        Data(E).Radiance = Radiance(:);
                        Data(E).Row = Row(:);
                        Data(E).Swath = Swath(:);
                    end
                end
            end
        cd J:\New_Retrieval_sp_mats
        savename=[satellite,'_',retrieval,'_',year,month,day];
        save(savename, 'Data')
        %clear Data 
        pack
        toc
        end
    end
end

%..........................................................................
%   
%  Calulating AMF used by nasa folks....
%
%..........................................................................

%{
%Calculate ColumnAmountNO2Polluted (Stratospheric NO2 Column)

for a=1:length(Data)
   b=find(Data(a).ColumnAmountNO2Polluted < 0);
   Data(a).ColumnAmountNO2Polluted(b)=0;
   Data(a).ColumnAmountNO2Unpolluted = Data(a).ColumnAmountNO2Initial-(Data(a).ColumnAmountNO2Polluted.*Data(a).AMFPolluted./Data(a).AMFUnpolluted);


% Calculate AMF applied for Standard Product NO2 Columns (AMFap)

    Data(a).AMFap = (Data(a).SlantColumnAmountNO2 - Data(a).ColumnAmountNO2Unpolluted .* Data(a).AMFUnpolluted) ./ Data(a).ColumnAmountNO2Polluted;
    c=find(Data(a).AMFap<0);
    Data(a).AMFap(c)=Data(a).AMFInitial(c);
end


%SLANT = Data(a).ColumnAmountNO2Polluted .* Data(a).AMFPolluted + Data(a).ColumnAmountNO2Unpolluted .* Data(a).AMFUnpolluted;
%}

%{
..........................................................................
OTHER STUFF

stride = [1 1 1];
blocksize = [1 1 1];
offset = [0,0];
slabsize = [10,10];
fileID = H5F.open(filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');

datasetID = H5D.open(fileID, hinfo.GroupHierarchy.Groups(1).Groups(2).Groups(1).Groups(2).Datasets(2).Name);
dataspaceID = H5D.get_space(datasetID);
H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', [1,1], stride, [10,10], blocksize);
memspaceID = H5S.create_simple(length(slabsize), slabsize, slabsize);
data = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT');
%}

