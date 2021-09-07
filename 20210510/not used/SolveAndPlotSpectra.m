%read and process field data
function SolveAndPlotSpectra
d=dir('*.sig');
addpath('/home/nbair/code/ParBal/');
addpath(genpath('/home/snowhydro/jdozier/MATLAB/SPIReS2021jd'));
addpath('/home/snowhydro/nbair/CUES/CUESspectrometer')
for i=1:length(d)
    fname=d(i).name;
    [T,~,S]=readSVC(fname);
    S.matdateUTC=S.matdateUTC;
    [declin,~,solar_lon]=Ephemeris(S.matdateUTC);
    cosZ=sunang2(S.lat,S.lon,declin,solar_lon);

    R0=table(T.wavelength,zeros(height(T),1),'VariableNames',...
        {'wavelength','reflectance'});
    R0.Properties.VariableUnits={'nm',''};
    
    rfit=fit(T.wavelength,T.ReflectanceTarget,'linearinterp');
%     P = setPrescription('snow','wavelength',T.wavelength,...
%         'r0',R0,'cosZ',cosZ,'waveUnit','nm','LAPradius',3,'LAP','dust',...
%         'LAPfraction',0,'elevation',3e3,'lookup',true);
%     unknowns={'fSCA','radius','lapfraction','wetness'};
%     load('atmosWeightSVC.mat','Ftrans');
%     
%     [o,ostats,P] = SPIReS_inv(rfit,'snow',unknowns,P,...
%         'atmosphere',Ftrans,'method','normResiduals');


    P = setPrescription('snow','wavelength',T.wavelength,...
        'r0',R0,'cosZ',cosZ,'waveUnit','nm','LAPradius',3,'LAP','dust',...
        'LAPfraction',0,'elevation',3e3,'lookup',true);
    unknowns={'fSCA','lapfraction','wetness'};
    load('atmosWeightSVC.mat','Ftrans');
    t=Ftrans.GridVectors{1}<1000;
    %wavelength (nm),altit,cosZ
    Ftrans.Values(t,:,:)=0;
    [o,ostats,P] = SPIReS_inv(rfit,'snow',unknowns,P,...
        'atmosphere',Ftrans,'method','normResiduals');
    o.radius=P.snow.radius;
    
    snow_str=sprintf(...
        ['modeled:\nfshade=%0.2f\nradius=%4.0f ',...
        'um\ndust conc=%3.0f ppmw,radius=%3.1f (given) um\n' ...
        'WE=%3.0f mm (given)\n wetness=%0.2f'],...
        1-o.fSCA,o.radius,o.LAPfraction*1e6,P.snow.LAPradius,...
        P.snow.waterEquivalent,o.wetness);
    
    if i==1
        tit=sprintf('Daves wave unsmoothed %i',i);
    elseif i==2
        tit=sprintf('Daves wave smoothed');
    elseif i==3
        tit=sprintf('Upper Gold hill unsmoothed');
    elseif i==5
        tit=sprintf('Upper Gold hill smoothed');
    end
    
    model=SPIReS_fwd(P);
    figure;
    plot(T.wavelength,[T.ReflectanceTarget model],'LineWidth',2);
    
    legend('measured',snow_str,'Interpreter',...
        'tex');
    title(['MMSA ' tit ' ' datestr(S.matdateUTC-8/24)],...
        'Interpreter','none');
    text(0,0.2,fname,'units','normalized','Interpreter','none');
    ylabel('reflectance');
    xlabel('wavelength, nm');
    ylim([0 1]);
end
