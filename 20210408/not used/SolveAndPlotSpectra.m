%read and process field data
function SolveAndPlotSpectra
solutionMethod='lsqnonlin';
d=dir('*.sig');
addpath('/home/nbair/code/ParBal/');

for i=1:length(d)
    fname=d(i).name;
    [T,hdr,S]=readSVC(fname);

%     t=T.RadianceTarget/max(T.RadianceTarget) < 0.002;
%     T.ReflectanceTarget(t)=0;
    S.matdateUTC=S.matdateUTC-1; % off by 1 day
    [declin,~,solar_lon]=Ephemeris(S.matdateUTC);
    cosZ=sunang(S.lat,S.lon,declin,solar_lon);
    if i~=4 && i~=5
            [o,ostats,P]=invertSnowCloudSpectralRefl(...
            T.ReflectanceTarget,{'radius','dust','wetness'},'snow',...
            'wavelength',T.wavelength,...
            'cosZ',cosZ,'waveUnit','nm','solutionMethod',solutionMethod);
            o.waterEquivalent=Inf;
            snow_str=sprintf(...
                'modeled:\nradius=%4.0f um\ndust=%3.0f ppmw\nWE=%3.0f mm\n wetness=%0.2f',...
            o.radius,o.dust*1e6,o.waterEquivalent,o.wetness);
    elseif i==5
        %run with shade endmember
        R0=table(T.wavelength,zeros(height(T),1),'VariableNames',...
            {'wavelength','reflectance'});
            [o,ostats,P]=invertSnowCloudSpectralRefl(...
            T.ReflectanceTarget,{'fSCA','radius','dust','wetness'},'snow',...
            'wavelength',T.wavelength,'R0',R0,...
            'cosZ',cosZ,'waveUnit','nm','solutionMethod',solutionMethod);
            o.waterEquivalent=Inf;
            snow_str=sprintf(...
                'modeled:\nfshade=%0.2f\nradius=%4.0f um\ndust=%3.0f ppmw\nWE=%3.0f mm\n wetness=%0.2f',...
                1-o.fSCA,o.radius,o.dust*1e6,o.waterEquivalent,o.wetness);
    else
        continue
    end
if i==1
    tit='Dennys donwhill';
elseif i==2
    tit='Top of Santiago';
elseif i==3
    tit='Hemlock flats';
elseif i==5
    tit='Hemlock flats, undisturbed';
end
    model=SnowCloudSpectralRefl(P);
    figure;
    plot(T.wavelength,[T.ReflectanceTarget model.refl],'LineWidth',2);
    
   
    legend('measured',snow_str,'Interpreter',...
        'tex');
    title(['MMSA ' tit ' ' datestr(S.matdateUTC-8/24)],...
        'Interpreter','none');
    text(0,0.2,fname,'units','normalized','Interpreter','none');
    ylabel('reflectance');
    xlabel('wavelength, um');
   
end
