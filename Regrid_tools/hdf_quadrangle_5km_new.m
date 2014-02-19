%%hdf_quadrangle_5km_new
%%arr 11/19/2009
%%uses functions exchange_coord,round01,calcline,clip 
%JLL 2-14-2014: round01 does not seem to be used; the built-in matlab
%function round seems to be in place


%maxx=501; maxy=1201;

Time=zeros(maxx,maxy);
ViewingZenithAngle=zeros(maxx,maxy);
SolarZenithAngle=zeros(maxx,maxy);
ViewingAzimuthAngle=zeros(maxx,maxy);
SolarAzimuthAngle=zeros(maxx,maxy);
CloudFraction=zeros(maxx,maxy);
CloudRadianceFraction=zeros(maxx,maxy);
ColumnAmountNO2=zeros(maxx,maxy);
SlantColumnAmountNO2=zeros(maxx,maxy);
TerrainHeight=zeros(maxx,maxy);
TerrainPressure=zeros(maxx,maxy);
TerrainReflectivity=zeros(maxx,maxy);
vcdQualityFlags=zeros(maxx,maxy);
CloudPressure=zeros(maxx,maxy);
RelativeAzimuthAngle=zeros(maxx,maxy);
Latitude=zeros(maxx,maxy);
Longitude=zeros(maxx,maxy);
Pixel=zeros(maxx,maxy);
ColumnAmountNO2Trop=zeros(maxx,maxy);
Areaweight=zeros(maxx,maxy);
Count=zeros(maxx,maxy);
Count2=0;
GLOBETerpres=zeros(maxx,maxy);
MODISAlbedo=zeros(maxx,maxy);
BEHRAMFTrop=zeros(maxx,maxy);
BEHRColumnAmountNO2Trop=zeros(maxx,maxy);
MODISCloud=zeros(maxx,maxy);
Row=zeros(maxx,maxy);
Swath=zeros(maxx,maxy);
Area=zeros(maxx,maxy);
AMFTrop=zeros(maxx,maxy);
AMFStrat=zeros(maxx,maxy);
ColumnAmountNO2Strat=zeros(maxx,maxy);
ColumnAmountNO2Initial=zeros(maxx,maxy);
XTrackQualityFlags=zeros(maxx,maxy);

for x=1:1:Dimensions(1)*Dimensions(2);
    num=x;
    y=1;

    pixelarea=(m_lldist([Lon1(x,y)-180 Lon2(x,y)-180],[Lat1(x,y) Lat2(x,y)])/1000)*13;

    x1=round(lCoordLat(x,y,1)); x1=clip(x1,minx,maxx);
    y1=round(lCoordLon(x,y,1)); y1=clip(y1,miny,maxy);
    x2=round(lCoordLat(x,y,2)); x2=clip(x2,minx,maxx);
    y2=round(lCoordLon(x,y,2)); y2=clip(y2,miny,maxy);
    x3=round(lCoordLat(x,y,3)); x3=clip(x3,minx,maxx);
    y3=round(lCoordLon(x,y,3)); y3=clip(y3,miny,maxy);
    x4=round(lCoordLat(x,y,4)); x4=clip(x4,minx,maxx);
    y4=round(lCoordLon(x,y,4)); y4=clip(y4,miny,maxy);

    
    if y2<y1; [x1,y1,x2,y2]=exchange_coord(x1,y1,x2,y2); end
    if y3<y1; [x1,y1,x3,y3]=exchange_coord(x1,y1,x3,y3); end
    if y4<y1; [x1,y1,x4,y4]=exchange_coord(x1,y1,x4,y4); end
    if y2>y3; [x2,y2,x3,y3]=exchange_coord(x2,y2,x3,y3); end
    if y4>y3; [x4,y4,x3,y3]=exchange_coord(x4,y4,x3,y3); end
    if x4>x2; [x4,y4,x2,y2]=exchange_coord(x4,y4,x2,y2); end
       
    
    
    Time_i=Data(d).Time;
    ViewingZenithAngle_i=Data(d).ViewingZenithAngle;
    SolarZenithAngle_i=Data(d).SolarZenithAngle;
    ViewingAzimuthAngle_i=Data(d).ViewingAzimuthAngle;
    SolarAzimuthAngle_i=Data(d).SolarAzimuthAngle;
    CloudFraction_i=Data(d).CloudFraction;
    CloudRadianceFraction_i=Data(d).CloudRadianceFraction;
    ColumnAmountNO2_i=Data(d).ColumnAmountNO2;
    ColumnAmountNO2Initial_i=Data(d).ColumnAmountNO2Initial;
    SlantColumnAmountNO2_i=Data(d).SlantColumnAmountNO2;
    TerrainHeight_i=Data(d).TerrainHeight;
    TerrainPressure_i=Data(d).TerrainPressure;
    TerrainReflectivity_i=Data(d).TerrainReflectivity;
    vcdQualityFlags_i=Data(d).vcdQualityFlags;
    CloudPressure_i=Data(d).CloudPressure;
    RelativeAzimuthAngle_i=Data(d).RelativeAzimuthAngle;
    Latitude_i=Data(d).Latitude;
    Longitude_i=Data(d).Longitude;
    Pixel_i=repmat(1:length(Data(d).Longitude),60,1)';
    ColumnAmountNO2Trop_i=Data(d).ColumnAmountNO2Trop;
    GLOBETerpres_i=Data(d).GLOBETerpres;
    MODISAlbedo_i=Data(d).MODISAlbedo;
    BEHRAMFTrop_i=Data(d).BEHRAMFTrop;
    BEHRColumnAmountNO2Trop_i=Data(d).BEHRColumnAmountNO2Trop;
    MODISCloud_i=Data(d).MODISCloud;
    Row_i=Data(d).Row;
    Swath_i=Data(d).Swath;
    AMFTrop_i=Data(d).AMFTrop;
    AMFStrat_i=Data(d).AMFStrat;
    ColumnAmountNO2Strat_i=Data(d).ColumnAmountNO2Strat;
    ColumnAmountNO2Initial_i=Data(d).ColumnAmountNO2Initial;
    XTrackQualityFlags_i=Data(d).XTrackQualityFlags;
    
    
    Time_val=Time_i(x);
    ViewingZenithAngle_val=ViewingZenithAngle_i(x);
    SolarZenithAngle_val=SolarZenithAngle_i(x);
    ViewingAzimuthAngle_val=ViewingAzimuthAngle_i(x);
    SolarAzimuthAngle_val=SolarAzimuthAngle_i(x);
    CloudFraction_val=CloudFraction_i(x);
    CloudRadianceFraction_val=CloudRadianceFraction_i(x);
    ColumnAmountNO2_val=ColumnAmountNO2_i(x);
    SlantColumnAmountNO2_val=SlantColumnAmountNO2_i(x);
    TerrainHeight_val=TerrainHeight_i(x);
    TerrainPressure_val=TerrainPressure_i(x);
    TerrainReflectivity_val=TerrainReflectivity_i(x);
    vcdQualityFlags_val=vcdQualityFlags_i(x);
    CloudPressure_val=CloudPressure_i(x);
    RelativeAzimuthAngle_val=RelativeAzimuthAngle_i(x);
    Latitude_val=Latitude_i(x);
    Longitude_val=Longitude_i(x);
    Pixel_val=Pixel_i(x);
    ColumnAmountNO2Trop_val=ColumnAmountNO2Trop_i(x);
    GLOBETerpres_val=GLOBETerpres_i(x);
    MODISAlbedo_val=MODISAlbedo_i(x);
    BEHRAMFTrop_val=BEHRAMFTrop_i(x);
    BEHRColumnAmountNO2Trop_val=BEHRColumnAmountNO2Trop_i(x);
    MODISCloud_val=MODISCloud_i(x);
    Row_val=Row_i(x);
    Swath_val=Swath_i(x);
    AMFTrop_val=AMFTrop_i(x);
    AMFStrat_val=AMFStrat_i(x);
    ColumnAmountNO2Strat_val=ColumnAmountNO2Strat_i(x);
    ColumnAmountNO2Initial_val=ColumnAmountNO2Initial_i(x);
    XTrackQualityFlags_val=XTrackQualityFlags_i(x);
    
    
    dim=[maxx maxy];
    bottom=y1+1;%y1 
    top=y3;%y3-1
    if (bottom<dim(2)) && (top>=1);
        bottom=clip(bottom,1,dim(2));
        top=clip(top,1,dim(2));
    end

    %if bottom<=0 
    %elseif top<=0;
    %    continue
    %else

    for y=bottom:1:top;
        if (x2>=calcline(y2,x1,y1,x3,y3)) && (x4<=calcline(y4,x1,y1,x3,y3));
            if y<y4; 
                left=calcline(y,x1,y1,x4,y4);
            else 
                left=calcline(y,x4,y4,x3,y3);
            end
            if y<y2; 
                right=calcline(y,x1,y1,x2,y2);
            else 
                right=calcline(y,x2,y2,x3,y3);
            end
        else 
            left=calcline(y,x1,y1,x3,y3);
            if y2>y4; 
                [x4,y4,x2,y2]=exchange_coord(x4,y4,x2,y2);
            end
            if y<y2; 
                right=calcline(y,x1,y1,x2,y2);
            elseif y<y4; 
                right=calcline(y,x2,y2,x4,y4);
            else 
                right=calcline(y,x4,y4,x3,y3);
            end
            if (x2<=calcline(y2,x1,y1,x3,y3)) && (x4<=calcline(y4,x1,y1,x3,y3));
                placeholder=left; 
                left=right;
                right=placeholder;
            end
        end
        right=right-1;
        if left<=0;
        elseif (dim(1)>=right>=1) && (1<=left<dim(1));
            clip(left,1,dim(1));
            clip(right,1,dim(1));
            Count2=Count2+1;
            for x=left+1:right+1;
            %for x=left:right;
                if Time(x,y)~=0 && isnan(ColumnAmountNO2Trop(x,y))==0;
                    
                    Time(x,y)=mean([Time(x,y);Time_val]);
                    ViewingZenithAngle(x,y)=mean([ViewingZenithAngle(x,y),ViewingZenithAngle_val]);
                    SolarZenithAngle(x,y)=mean([SolarZenithAngle(x,y),SolarZenithAngle_val]);
                    ViewingAzimuthAngle(x,y)=mean([ViewingAzimuthAngle(x,y),ViewingAzimuthAngle_val]);
                    SolarAzimuthAngle(x,y)=mean([SolarAzimuthAngle(x,y),SolarAzimuthAngle_val]);
                    CloudFraction(x,y)=mean([CloudFraction(x,y),CloudFraction_val]);
                    CloudRadianceFraction(x,y)=mean([CloudRadianceFraction(x,y),CloudRadianceFraction_val]);
                    ColumnAmountNO2(x,y)=mean([ColumnAmountNO2(x,y),ColumnAmountNO2_val]);
                    SlantColumnAmountNO2(x,y)=mean([SlantColumnAmountNO2(x,y),SlantColumnAmountNO2_val]);
                    TerrainHeight(x,y)=mean([TerrainHeight(x,y),TerrainHeight_val]);
                    TerrainPressure(x,y)=mean([TerrainPressure(x,y),TerrainPressure_val]);
                    TerrainReflectivity(x,y)=mean([TerrainReflectivity(x,y),TerrainReflectivity_val]);
                    vcdQualityFlags(x,y)=mean([vcdQualityFlags(x,y),vcdQualityFlags_val]);
                    CloudPressure(x,y)=mean([CloudPressure(x,y),CloudPressure_val]);
                    RelativeAzimuthAngle(x,y)=mean([RelativeAzimuthAngle(x,y),RelativeAzimuthAngle_val]);
                    Latitude(x,y)=mean([Latitude(x,y),Latitude_val]);
                    Longitude(x,y)=mean([Longitude(x,y),Longitude_val]);
                    Pixel(x,y)=mean([Pixel(x,y),Pixel_val]);
                    Area(x,y)=mean([Area(x,y);pixelarea]);
                    Areaweight(x,y)=2/Area(x,y);
                    ColumnAmountNO2Trop(x,y)=mean([ColumnAmountNO2Trop(x,y),ColumnAmountNO2Trop_val]);
                    GLOBETerpres(x,y)=mean([GLOBETerpres(x,y),GLOBETerpres_val]);
                    MODISAlbedo(x,y)=mean([MODISAlbedo(x,y),MODISAlbedo_val]);
                    BEHRAMFTrop(x,y)=mean([BEHRAMFTrop(x,y),BEHRAMFTrop_val]);
                    BEHRColumnAmountNO2Trop(x,y)=mean([BEHRColumnAmountNO2Trop(x,y),BEHRColumnAmountNO2Trop_val]);
                    MODISCloud(x,y)=mean([MODISCloud(x,y),MODISCloud_val]);
                    Row(x,y)=mean([Row(x,y),Row_val]);
                    Swath(x,y)=mean([Swath(x,y),Swath_val]);
                    Count(x,y)=Count(x,y)+1;
                    AMFTrop(x,y)=mean([AMFTrop(x,y),AMFTrop_val]);
                    AMFStrat(x,y)=mean([AMFStrat(x,y),AMFStrat_val]);
                    ColumnAmountNO2Strat(x,y)=mean([ColumnAmountNO2Strat(x,y),ColumnAmountNO2Strat_val]);
                    ColumnAmountNO2Initial(x,y)=mean([ColumnAmountNO2Initial(x,y),ColumnAmountNO2Initial_val]);
                    XTrackQualityFlags(x,y)=mean([XTrackQualityFlags(x,y),XTrackQualityFlags_val]);
                else               
                    Time(x,y)=Time_val;
                    ViewingZenithAngle(x,y)=ViewingZenithAngle_val;
                    SolarZenithAngle(x,y)=SolarZenithAngle_val;
                    ViewingAzimuthAngle(x,y)=ViewingAzimuthAngle_val;
                    SolarAzimuthAngle(x,y)=SolarAzimuthAngle_val;
                    CloudFraction(x,y)=CloudFraction_val;
                    CloudRadianceFraction(x,y)=CloudRadianceFraction_val;
                    ColumnAmountNO2(x,y)=ColumnAmountNO2_val;
                    SlantColumnAmountNO2(x,y)=SlantColumnAmountNO2_val;
                    TerrainHeight(x,y)=TerrainHeight_val;
                    TerrainPressure(x,y)=TerrainPressure_val;
                    TerrainReflectivity(x,y)=TerrainReflectivity_val;
                    vcdQualityFlags(x,y)=vcdQualityFlags_val;
                    CloudPressure(x,y)=CloudPressure_val;
                    RelativeAzimuthAngle(x,y)=RelativeAzimuthAngle_val;
                    Latitude(x,y)=Latitude_val;
                    Longitude(x,y)=Longitude_val;
                    Pixel(x,y)=Pixel_val;
                    Area(x,y)=pixelarea;
                    Areaweight(x,y)=1/Area(x,y);
                    ColumnAmountNO2Trop(x,y)=ColumnAmountNO2Trop_val;
                    GLOBETerpres(x,y)=GLOBETerpres_val;
                    MODISAlbedo(x,y)=MODISAlbedo_val;
                    BEHRAMFTrop(x,y)=BEHRAMFTrop_val;
                    BEHRColumnAmountNO2Trop(x,y)=BEHRColumnAmountNO2Trop_val;
                    MODISCloud(x,y)=MODISCloud_val;
                    Row(x,y)=Row_val;
                    Swath(x,y)=Swath_val;
                    Count(x,y)=Count(x,y)+1;
                    AMFTrop(x,y)=AMFTrop_val;
                    AMFStrat(x,y)=AMFStrat_val;
                    ColumnAmountNO2Strat(x,y)=ColumnAmountNO2Strat_val;
                    ColumnAmountNO2Initial(x,y)=ColumnAmountNO2Initial_val;
                    XTrackQualityFlags(x,y)=XTrackQualityFlags_val;
                end
            end
        end
    end
end 
