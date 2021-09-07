function Tbl = runner(folder,varargin)
% Tbl = runner(folder [,doPlot]
p = inputParser;
addRequired(p,'folder',@ischar);
addOptional(p,'doPlot',false,@(x) isnumeric(x) || islogical(x));
parse(p,folder,varargin{:});
doPlot = p.Results.doPlot;
unkf = {'radius','wetness','fSCA','lapfraction'};
unkn = {'radius','wetness','lapfraction'};
f = dir(fullfile(folder,'2*.mat'));
load('atmosWeightSVC.mat','Ftrans')
Tbl = table;
modelType = categorical({'normR shade', 'norm residuals','least squares','spec ang'});
nMod = length(modelType);
for k=1:length(f)
    load(fullfile(f(k).folder,f(k).name),'S','T')
    if ~isempty(S.lat)
        [declin,~,omega] = EarthEphemeris(S.datetime);
        cosZ = sunang(S.lat,S.lon,declin,omega);
    else
        cosZ = 2/3;
    end
    P0 = setPrescription('snow','wavelength',T.wavelength,'waveunit','nm',...
        'wetness',.02,'lapfraction',[1e-8 1e-10],'LAP',{'dust','soot'},'lapradius',[5 .25],'lookup',true,'elevation',2940,'cosZ',cosZ);
    Rfun = griddedInterpolant(T.wavelength,T.ReflectanceTarget,'makima');
    o = cell(nMod,1);
    s = cell(size(o));
    p = cell(size(o));
    useModel = false(nMod,1);
    for r=1:nMod
        switch r
            case 1
                % fmincon, solve for fSCA too
                unk = unkf;
                method = 'norm';
            case 2
                % fmincon, assume fSCA = 1
                method = 'norm';
                unk = unkn;
            case 3
                % least squares, assume fSCA = 1
                unk = unkn;
                method = 'lsq';
            case 4
                % fmincon, spectral angle
                unk = unkn;
                method = 'spec';
        end
        [o{r},s{r},p{r}] = SPIReS_inv(Rfun,'snow',unk,P0,'atmos',Ftrans,'method',method);
        
        if s{r}.exitflag>0
            useModel(r) = true;
            if r==4
                s{r}.normResiduals = s{r}.sineAngle;
            end
        end
    end
    
    % use all in table except those with bad exitflag
    [~,tn,~] = fileparts(f(k).name);
    for r=1:nMod
        if useModel(r)
            thisTbl = table(categorical({tn}),S.datetime,categorical({o{r}.solver}),...
                modelType(r),p{r}.snow.radius,p{r}.snow.wetness,p{r}.snow.LAPfraction,p{r}.snow.fSCA,...
                s{r}.normResiduals,s{r}.slope,s{r}.rmse,s{r}.goodness,...
                'VariableNames',{'file','datetime','solver','model','radius','wetness',...
                'LAPfraction','fSCA','normRsinA','slope','RMSE','goodness'});
            Tbl = [Tbl; thisTbl]; %#ok<AGROW>
        end
    end
    if doPlot
        figure
        plot(T.wavelength,T.ReflectanceTarget,'k--','linewidth',1.5)
        hold on;
        lineColor = {'b','k','r','g'};
        for r=1:nMod
            if useModel(r)
                plot(T.wavelength,SPIReS_fwd(p{r}),'linewidth',1,'color',lineColor{r})
            end
        end
        title(tn,'Interpreter','none')
        xlabel('wavelength (nm)')
        legText = {'measure'};
        for r=1:nMod
            if useModel(r)
                legText = cat(2,legText,char(modelType(r)));
            end
        end
        legend(legText);
        axis padded
        saveas(gcf,[tn '.png'])
    end
end
end