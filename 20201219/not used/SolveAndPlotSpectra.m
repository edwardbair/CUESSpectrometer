%read and process field data
function SolveAndPlotSpectra
solutionMethod='lsqnonlin';
d=dir('*.sig');

for i=1:length(d)
    fname=d(i).name;
    [T,hdr,S]=readSVC(fname);
%     t=T.RadianceTarget<75;
%     T.ReflectanceTarget(t)=0;
    [declin,~,solar_lon]=Ephemeris(S.matdateUTC);
    cosZ=sunang(S.lat,S.lon,declin,solar_lon);
    if i==1 || i==3 %snow only
        [o,ostats,P]=invertSnowCloudSpectralRefl(...
            T.ReflectanceTarget,{'radius','dust'},'snow',...
            'wavelength',T.wavelength,...
            'cosZ',cosZ,'waveUnit','nm','solutionMethod',solutionMethod);
            o.waterEquivalent=Inf;
            tit='thick snow';
            o.wetness=0;
    elseif i==2 % dark target
      
        R0=readSVC(d(4).name);
        R0.Properties.VariableNames{4}='reflectance';
        R0.reflectance(R0.reflectance<0)=0;
        [o,ostats,P]=invertSnowCloudSpectralRefl(...
            T.ReflectanceTarget,{'radius','waterEquivalent'},'snow',...
            'wavelength',T.wavelength,'R0',R0,...
            'cosZ',cosZ,'waveUnit','nm','solutionMethod',solutionMethod);
        o.dust=0;
            tit='thin snow, gray target';
            o.wetness=0;
    elseif i==4
        continue;
    elseif i==5 || i==6
            [o,ostats,P]=invertSnowCloudSpectralRefl(...
            T.ReflectanceTarget,{'radius','dust','wetness'},'snow',...
            'wavelength',T.wavelength,...
            'cosZ',cosZ,'waveUnit','nm','solutionMethod',solutionMethod);
            o.waterEquivalent=Inf;
            tit='dirty wet snow';
    end
    model=SnowCloudSpectralRefl(P);
    figure;
    
    plot(T.wavelength,[T.ReflectanceTarget model.refl],'LineWidth',2);
    legend('measured',...
        sprintf('modeled:\nradius=%4.0f um\ndust=%3.0f ppmw\nWE=%3.0f mm\n wetness=%0.2f',...
        o.radius,o.dust*1e6,o.waterEquivalent,o.wetness),'Interpreter',...
        'tex');
    title(['CUES ' tit ' ' datestr(S.matdateUTC-8/24)],...
        'Interpreter','none');
    text(0,0.2,fname,'units','normalized','Interpreter','none');
    ylabel('reflectance');
    xlabel('wavelength, um');
   
end
