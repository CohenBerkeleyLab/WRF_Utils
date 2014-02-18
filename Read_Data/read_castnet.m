%read_castnet
%arr 12/14/2011

addpath G:\CASTNET

years={'2005';'2006';'2007';'2008';'2009';'2010';'2011'};

nlines=1;
for y=1:length(years)
    year=years{y};
    
    filename=['metdata_',year,'.csv'];

    howbig=1E6;  % guess how many rows of data and grab the space first

    cnSiteID=cell(1,howbig); 
    cnDay=zeros(1,howbig); cnMonth=zeros(1,howbig);
    cnYear=zeros(1,howbig); cnHour=zeros(1,howbig);
    cnTemperature=zeros(1,howbig); cnTemperatureF=cell(1,howbig);
    cnRH=zeros(1,howbig); cnRHF=cell(1,howbig);
    cnSolarRad=zeros(1,howbig); cnSolarRadF=cell(1,howbig);
    cnOzone=zeros(1,howbig); cnOzoneF=cell(1,howbig);
    cnPrecip=zeros(1,howbig); cnPrecipF=cell(1,howbig);
    cnWindSpeed=zeros(1,howbig); cnWindSpeedF=cell(1,howbig);
    
    fid=fopen(filename, 'r');
    line=fgets(fid);
    line=fgets(fid);
    while line~=-1
    text = textscan(line, '%q', 'delimiter', ',');
    text_char = char(text{1});

    thisSiteID=text_char(1,:);
       if isempty(thisSiteID);
          thisSiteID=NaN;
       end
    cnSiteID{nlines}=thisSiteID;

    thisDate=text_char(2,:);
       if isempty(thisDate);
          thisDate=NaN;
       end
    cnDay(nlines)=str2double(thisDate(9:10));
    cnMonth(nlines)=str2double(thisDate(6:7));
    cnYear(nlines)=str2double(thisDate(1:4));
    cnHour(nlines)=str2double(thisDate(12:13));
        
    thisTemperature=str2double(text_char(3,:));
       if isempty(thisTemperature);
          thisTemperature=NaN;
       end
    cnTemperature(nlines)=thisTemperature;
        
    thisTemperatureF=text_char(4,:);
       if isempty(thisTemperatureF);
          thisTemperatureF=NaN;
       end
    cnTemperatureF{nlines}=thisTemperatureF;
        
    thisRH=str2double(text_char(7));
       if isempty(thisRH);
          thisRH=NaN;
       end
    cnRH(nlines)=thisRH;

    thisRHF=text_char(8);
       if isempty(thisRHF);
          thisRHF=NaN;
       end
    cnRHF{nlines}=thisRHF;

    thisSolarRad=str2double(text_char(9));
       if isempty(thisSolarRad);
          thisSolarRAd=NaN;
       end
    cnSolarRad(nlines)=thisSolarRad;

    thisSolarRadF=text_char(10);
       if isempty(thisSolarRadF);
          thisSolarRadF=NaN;
       end
    cnSolarRadF{nlines}=thisSolarRadF;

    thisOzone=str2double(text_char(11));
       if isempty(thisOzone);
          thisOzone=NaN;
       end
    cnOzone(nlines)=thisOzone;
       
    thisOzoneF=text_char(12);
       if isempty(thisOzoneF);
          thisOzoneF=NaN;
       end
    cnOzoneF{nlines}=thisOzoneF;

    thisPrecip=str2double(text_char(13));
       if isempty(thisPrecip);
          thisPrecip=NaN;
       end
    cnPrecip(nlines)=thisPrecip;

    thisPrecipF=text_char(14);
       if isempty(thisPrecipF);
          thisPrecipF=NaN;
       end
    cnPrecipF{nlines}=thisPrecipF;
        
    thisWindSpeed=str2double(text_char(15));
        if isempty(thisWindSpeed);
           thisWindSpeed=NaN;
        end
    cnWindSpeed(nlines)=thisWindSpeed;

    thisWindSpeedF=text_char(16);
        if isempty(thisWindSpeedF);
           thisWindSpeedF=NaN;
        end
    cnWindSpeedF{nlines}=thisWindSpeedF;
    
        line=fgets(fid);
        nlines=nlines+1;
       if mod(nlines,2000)==0
          disp([num2str(nlines),' lines read so far'])
       end
    end
fclose(fid); 


cnSiteID=cnSiteID(1:nlines); cnSiteID=deblank(cnSiteID);
cnDay=cnDay(1:nlines);
cnMonth=cnMonth(1:nlines); 
cnYear=cnYear(1:nlines); 
cnHour=cnHour(1:nlines); 
cnTemperature=cnTemperature(1:nlines); 
cnTemperatureF=cnTemperatureF(1:nlines); cnTemperatureF=deblank(cnTemperatureF);
cnRH=cnRH(1:nlines);  
cnRHF=cnRHF(1:nlines); cnRHF=deblank(cnRHF);
cnSolarRad=cnSolarRad(1:nlines); 
cnSolarRadF=cnSolarRadF(1:nlines); cnSolarRadF=deblank(cnSolarRadF);
cnOzone=cnOzone(1:nlines); 
cnOzoneF=cnOzoneF(1:nlines); cnOzoneF=deblank(cnOzoneF);
cnPrecip=cnPrecip(1:nlines); 
cnPrecipF=cnPrecipF(1:nlines); cnPrecipF=deblank(cnPrecipF);
cnWindSpeed=cnWindSpeed(1:nlines); 
cnWindSpeedF=cnWindSpeedF(1:nlines); cnWindSpeedF=deblank(cnPrecipF);
  
cnWeekday=zeros(1,length(nlines));
for i=1:length(cnDay);
    day=num2str(cnDay(i));
    month=num2str(cnMonth(i));
    date=[month,'/',day,'/',year];
    a=datenum(date);
    [n, s] = weekday(a);
    cnWeekday(i)=n;
end

cd('G:/CASTNET/')
savename=['CASTNET',year,'data'];

save(savename, 'cnSiteID', 'cnDay', 'cnMonth', 'cnYear', 'cnHour', 'cnTemperature', 'cnTemperatureF', 'cnRH', 'cnRHF', 'cnSolarRad', 'cnSolarRadF', 'cnOzone', 'cnOzoneF', 'cnPrecip', 'cnPrecipF', 'cnWindSpeed', 'cnWindSpeedF');

end