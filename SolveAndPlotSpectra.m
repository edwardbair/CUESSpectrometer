function T=SolveAndPlotSpectra(loc)
%fit SVC spectra
%input: loc, directory to process

d=dir(fullfile(loc,'*.sig'));
addpath('/home/nbair/code/ParBal/');
addpath('/home/snowhydro/nbair/CUES/CUESspectrometer');

for i=1:length(d)
    fname=fullfile(d(i).folder,d(i).name);
    [T,~,S]=readSVC(fname);
    S.matdateUTC=S.matdateUTC;
    [declin,~,solar_lon]=Ephemeris(S.matdateUTC);
    cosZ=sunang2(S.lat,S.lon,declin,solar_lon);
    
    
    R0=table(T.wavelength,zeros(height(T),1),'VariableNames',...
        {'wavelength','reflectance'});
    R0.Properties.VariableUnits={'nm',''};
    
    rfit=fit(T.wavelength,T.ReflectanceTarget,'linearinterp');

         P = setPrescription('snow','wavelength',T.wavelength,...
        'r0',R0,'cosZ',cosZ,'waveUnit','nm','LAPradius',3,'LAP','SanJuanDust',...
        'LAPfraction',0,'elevation',3e3,'lookup',false);   
    
%     P = setPrescription('snow','wavelength',T.wavelength,...
%         'r0',R0,'cosZ',cosZ,'waveUnit','nm','LAPradius',3,'LAP','dust',...
%         'LAPfraction',0,'elevation',3e3,'lookup',true);
% RI=load('pumiceRefractiveIndex.mat');
% P = setPrescription('snow','wavelength',T.wavelength,...
%         'r0',R0,'cosZ',cosZ,'waveUnit','nm','LAPradius',3,...
%         'LAP',RI.VolcanicPumiceVolz73,...
%         'LAPfraction',0,'elevation',3e3,'lookup',false);

    unknowns={'fSCA','radius','lapfraction'};
    
    [o,ostats,P] = SPIReS_inv(rfit,'snow',unknowns,P,'method','lsq','useParallel',true);
    
    snow_str=sprintf(...
        ['modeled (rmse=%0.3f)\nfshade=%0.2f\nradius=%4.0f ',...
        'um\ndust conc=%3.0f ppmw,radius=%3.1f (given) um\n' ...
        ],...
        ostats.rmse,1-o.fSCA,o.radius,o.LAPfraction*1e6,...
        P.snow.LAPradius);
    
    model=SPIReS_fwd(P);
    figure;
    plot(T.wavelength,[T.ReflectanceTarget model],'LineWidth',2);
    
    legend('measured',snow_str,'Interpreter',...
        'tex');
    title([S.comment ' ' datestr(S.matdateUTC-8/24)],...
        'Interpreter','none');
    %     text(0,0.2,fname,'units','normalized','Interpreter','none');
    ylabel('reflectance');
    xlabel('wavelength, nm');
    ylim([0 1]);
end
