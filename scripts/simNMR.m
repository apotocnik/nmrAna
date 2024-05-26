%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT_NaLSX_1                         %
% Fitting Na12/Na-LSX spectra         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function simNMR(handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
COLOR = {[1,0,1],[1,0,0],[0,1,0],[0,0,1],[0.75,0.5,0]};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
if 1

    % Units
    e0 = 1.60217733e-19;
    h = 6.6260755e-34;
    tmpREF = [];

    %PATHname = 'F:\My Science\';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Na12/Na-LSX, powder 23Na at 400 MHz, Tdep
    A = []; B = []; freqL = []; S = []; f = [];
    param = get(handles.edtAnalyse,'String');
    eval(sprintf('%s\n',param{:}));
if SIM==1
    TEM = evalin('base','TEM');    

    % other parameters
    TITLE = ['Na_{12}/Na-LSX, ^{23}Na, powder, 400 MHz, ',num2str(TEM),'K'];
    PLOTtype = 1; Dy = 1; INDf = 1; Fsize = 8; normX = 1e3;
    SAVE = 1;
    saveNAME = ['Na12_test_',num2str(TEM),'K_M',num2str(MODEL),'_LB',num2str(LBtype),'_'];

    % experiment
    if 0 %peter script
        expF = 1;
        LB = 2e3;
        %expSPC = load([PATHname,'Na0_3_',num2str(LB/1e3),'_',num2str(TEM),'K_SMn.txt'],'-ascii');
        expSPC = load([PATHname,'Na0_2_3_',num2str(TEM),'K_SMn.txt'],'-ascii');
        
        %IND = find(expSPC(:,1)<=-0.1|expSPC(:,1)>=0.1);
        %expSPC = expSPC(IND,:);
        expSPC(:,1) = expSPC(:,1)*1e6;%-3.5e3;
        expSPC(:,2) = expF*expSPC(:,2)/max(expSPC(:,2));
        %expSPC(:,2) = expF*expSPC(:,2);
    else %mine
        expF = 1;
        LB = 2e3;
        
        %i should make global funktion for tihs:
        mode = evalin('base','mode');
        PATHname = evalin('base','PATHname');
        FileName = get(handles.edtFN,'String');
[~, Name Ext] = fileparts(FileName);
    if evalin('base','exist(''LSpar'',''var'')')
        LSpar = evalin('base', 'LSpar');
    else
        LSpar = {};
        LSpar{1}  = 0;
    end        
    if evalin('base','exist(''LSname'',''var'')')
        LSname = evalin('base','LSname');
    else
        LSname = 'noname';
    end

     if LSpar{1} == 1
            ind = str2double(LSpar{3});
            if (ind==0)
             	cla(handles.axeRes);
                cla(handles.axePhi);
                cla(handles.axeSHL);
                cla(handles.axeSPC);
                return;
            end
            splitName = regexp(Name,mode,'split');
            nameIN = ['lineShift\' char(splitName(1)) 'LS-' LSname '-' LSpar{3} '_'];
            if evalin('base','exist(''normal'',''var'')')
                normal = evalin('base','normal');
            end
    elseif strcmp(mode,'SW') 
        if evalin('base','exist(''SWname'',''var'')')
            SWname = evalin('base','SWname');
        else
            SWname = 'noname';
        end
        splitName = regexp(Name,mode,'split');
        nameIN = ['lineShift\' char(splitName(1)) 'SW-' SWname '_'];
    else
        splitName = regexp(Name,mode,'split');
        nameIN = ['lineShift\' char(splitName(1)) mode '_'];
     end
        
        expSPC = load([PATHname,nameIN,num2str(TEM),'K.txt'],'-ascii');
        REF = evalin('base','REF');
        
        expSPC(:,1) = expSPC(:,1)/1e6*REF*1e6;%-3.5e3;
        expSPC(:,2) = expF*expSPC(:,2)/max(expSPC(:,2)); 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MODEL
% allpar(1) = Vzz1
% allpar(2) = eta1
% allpar(3) = iso1
% allpar(4) = aksq1
% allpar(5) = akseta1
% allpar(6) = sig1
% allpar(7) = A1
% ...
% allpar(7*NS+1) = B

if MODEL==1
    %indpar = [7,14,21,28,35];
    indpar = [6,7,13,14,20,21,27,28,34,35];
    %indpar = [1,2,7,8,9,14];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Distribution of tensors

fi = [];

for i=1:length(freqL)

    % Powder (spherical ZCW)
    if DIStype{i}==1
        %N = 6044; g1 = 1427; g2 = 1891;
        N = 20000; g1 = 37; g2 = 199;
        %N = 50000; g1 = 1427; g2 = 1891;
        %N = 500000; g1 = 1427; g2 = 1891;
        UNI = (0:N-1)/N;
        phi = 2*pi*mod(UNI*g1,1);
        costheta = 2*mod(UNI*g2,1)-1;
        %costheta = mod(UNI*g2,1);
        theta = acos(costheta);
        sintheta = sin(theta);
    end

    fi{i} = [cos(2*phi);costheta.^2;costheta;sintheta];
    Rfi{i} = [cos(phi).*sin(theta);sin(phi).*sin(theta);cos(theta)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fitting function
if 1
    %frequency
    f = (-1:0.001:1)'*fR;

    %exp
    fEXP = expSPC(:,1)-tmpREF*freqL{1};
    fIND = find(fEXP>min(f)&fEXP<max(f));
    fEXP = fEXP(fIND);
    expSPC = expSPC(fIND,2);

    % fit
    if 1
        GAUSS = [];
        if FITtype==1
            if ismember(MODEL,[1:9])
                allinitpar = [];
                for i=1:NS
                    allinitpar =[allinitpar,Vzz(i),eta(i),iso(i),aksq(i),akseta(i),sig(i),A(i)];
                end
                allinitpar = [allinitpar,B];
            end
        end
        if FITtype==2
            allinitpar = allfitpar;
        end
        if FITtype==3
            allinitpar =load([PATHname,saveNAME,'par2.txt'],'savePAR2','-ascii');
        end
        if ismember(FITtype,[1,2,3])
            [fitpar,model] = simNMR_MODEL_1(f,fEXP,expSPC,freqL,fi,Rfi,tmpREF,Qmom,S,LBtype,MODEL,Niter,allinitpar,indpar,GAUSS,NS);
        end
        [sse,fitSPC,simSPC,simSPCall] = model(fitpar);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PLOT
if 1
    a4 = handles.axeRes;
    cla;
    hold(a4,'on');
    xlim(a4,'auto');ylim(a4,'auto');

    if PLOTtype==1
        allfitpar = allinitpar;
        allfitpar(indpar) = fitpar;
        varpar = 0.6*ones(size(allfitpar));
        varpar(indpar) = 0;

        plot(a4,fEXP/normX,expSPC,'-','linewidth',2,'Color',[1,0,0])
        %plot(a4,fEXP/normX,expSPC,'o','linewidth',2,'Color',[1,0,0],'MarkerSize',3)
        plot(a4,f/normX,simSPCall,'-','linewidth',2,'Color',[1,1,1]*0);
        for i=1:NS
            plot(a4,f/normX,simSPC{i},'-','linewidth',1,'Color',COLOR{i});
        end

        for i=1:NS
            if ismember(MODEL,[1:9])
                fit_Vzz = allfitpar(1+7*(i-1));
                fit_eta = allfitpar(2+7*(i-1));
                fit_iso = allfitpar(3+7*(i-1));
                fit_aksq = allfitpar(4+7*(i-1));
                fit_akseta = allfitpar(5+7*(i-1));
                fit_sig = allfitpar(6+7*(i-1));
                fit_A = allfitpar(7+7*(i-1));
                fit_freqQ = 3*e0*Qmom{i}*fit_Vzz/(h*(2*S{i}*(2*S{i}-1)));
                if i==1
                    fit_B = allfitpar(NS*7+1);
                end

                text(a4,0.02+0.2*(i-1),0.95,['{\itV_{ZZ}} =',num2str(round(1000*fit_Vzz/1e21)/1000),'x10^{21}V/m^2'],'FontSize',Fsize,'FontWeight','bold','Color',COLOR{i}*varpar(1),'Units','Norm')
                text(a4,0.02+0.2*(i-1),0.91,['{\it\eta} =',num2str(round(1000*fit_eta)/1000)],'FontSize',Fsize,'FontWeight','bold','Color',COLOR{i}*varpar(2),'Units','Norm')
                text(a4,0.02+0.2*(i-1),0.87,['{\it\delta}_{iso} =',num2str(round(fit_iso*1e6)),'ppm'],'FontSize',Fsize,'FontWeight','bold','Color',COLOR{i}*varpar(3),'Units','Norm')
                text(a4,0.02+0.2*(i-1),0.83,['{\it\delta}_{aniso} =',num2str(round(fit_aksq*1e6)),'ppm'],'FontSize',Fsize,'FontWeight','bold','Color',COLOR{i}*varpar(4),'Units','Norm')
                text(a4,0.02+0.2*(i-1),0.79,['{\it\eta}_{aniso} =',num2str(round(1000*fit_akseta)/1000)],'FontSize',Fsize,'FontWeight','bold','Color',COLOR{i}*varpar(5),'Units','Norm')
                text(a4,0.02+0.2*(i-1),0.75,['{\it\sigma} =',num2str(round(10*fit_sig/1e3)/10),'kHz'],'FontSize',Fsize,'FontWeight','bold','Color',COLOR{i}*varpar(7),'Units','Norm')
                text(a4,0.02+0.2*(i-1),0.71,['{\itA} =',num2str(round(1e6*fit_A)/1e6)],'FontSize',Fsize,'FontWeight','bold','Units','Norm')
                text(a4,0.02+0.2*(i-1),0.67,['{\it\nu}_Q =',num2str(round(1000*fit_freqQ/1e6)/1000),'MHz'],'FontSize',Fsize,'FontWeight','bold','Color',COLOR{i}*varpar(1),'Units','Norm')

                if i==1
                    text(a4,0.02+0.2*(i-1),0.40,['{\itB} =',num2str(round(1000*fit_B)/1000)],'FontSize',Fsize,'FontWeight','bold','Units','Norm')
                    text(a4,0.02+0.2*(i-1),0.36,['{\Delta}{^2} =',num2str(round(1000*sse)/1000)],'FontSize',Fsize,'FontWeight','bold','Units','Norm')
                end
            end
        end
        if normX==1e6
            xlabel(a4,'\nu-\nu_{ref} (MHz)'); ylabel(a4,'Intensity');
        elseif normX==1e3
            xlabel(a4,'\nu-\nu_{ref} (kHz)'); ylabel(a4,'Intensity');
        end
        title([TITLE,', MODEL = ',num2str(MODEL)]);
        %axis([min(f{INDf}/normX),max(f{INDf}/normX),-0.1*Dy,1.1*Dy+(i-1)*Dy]);
        %legend('exp','fit','partial',2)
        legend('exp','fit',1)
        %set(gca,'XDir','reverse')
        hold(a4,'off');
        %grid on;
        drawnow
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SAVE - remove?
if 0 %SAVE==1
    savePAR = [allfitpar,sse]';
    savePAR2 = allfitpar;

    if 1
        save([PATHname,saveNAME,'par.txt'],'savePAR','-ascii');
        save([PATHname,saveNAME,'par2.txt'],'savePAR2','-ascii');
    end
    if 1
        saveEXP = [(fEXP+tmpREF*freqL{1})/normX,expSPC];
        saveFIT = [(f+tmpREF*freqL{1})/normX,simSPCall];
        save([PATHname,saveNAME,'exp.txt'],'saveEXP','-ascii');
        save([PATHname,saveNAME,'fit.txt'],'saveFIT','-ascii');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    [-fit_freqQ/2*(1-fit_eta),-fit_freqQ/2*(1+fit_eta),fit_freqQ]/1e6
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
