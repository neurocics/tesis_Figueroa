% Time-frequiency calculations 
clear
% eeg analisis 
Sujetos = {
'AGBB_26121972' % 'RL WM GN'     % OK
'ASCS_31121976' % 'RL WM GN'  % sp1
'BBMA_19041977' % 'RL GN WM'  % sp1
'BCUV_02091982' % RL GN WM'  % sp1
'C_DG_15111995' % RL GN WM % sp1 GN
'CACB_27091965' % RL GN WM % sp1 
'CAMO_08111988' % 'RL GN WM' 
'CARM_14011974' % 'RL GN WM'
'CDPR_25031986' % 'RL GN WM'
'CEBT_30091989' % 'RL GN WM' % sp1 
'CJGV_24021957' % 'RL GN WM' % sp1 
'CMCO_17051975' % 'RL GN WM' 
'COGC_18022000' % 'RL GN WM' 
'EERV_01041955' % 'RL GN WM' 
'ELFN_09021961' % 'RL GN WM' 
'FABN_19111984' % 'RL GN WM' 
'FJAM_10012001' % 'RL GN WM' % sp1 
'FRCA_12081956' % 'RL GN WM' % sp1 
'GADM_08111983' % 'RL GN WM' % sp1 
'GAQJ_21061987' % 'RL GN WM' 
'IABF_16091997' % 'RL GN WM' 
'IABM_03061982' % 'RL GN WM' %RRR
'JAPV_03111995' % 'RL GN WM'
'JAZV_27081991' % 'RL GN WM'
'JCZD_03071963' % 'RL GN WM'% RRR % sp1 
'JDRP_06081956' % 'RL GN WM'  % sp1 
'JFFA_04101976' % 'RL GN WM'  % sp1 
'JFSO_05061967' % 'RL GN WM'  % sp1 
'JGRF_24091991' % 'RL GN WM'  % sp1 
'LAAB_29101982' % 'RL GN WM'  % sp1
'LDEV_31011988' % 'RL GN WM'  % sp1
'LMRG_21041962' % 'RL GN WM'  % sp1
'MAEM_17051957' % 'RL GN WM'  % sp1
'MFCM_21011984' % 'RL GN WM'  % sp1
'MGRN_22101979' % 'RL GN WM'  % sp1
'MLAD_09021957' % 'RL GN WM'  % sp1
'MLCB_11041981' % 'RL GN WM'  % sp1
'MPMM_22111977' % 'RL GN WM'  % sp1
'MTVV_27111974' % 'RL GN WM'  % sp1
'OAPC_01041962' % 'RL GN WM'  % sp1
'RABG_09091981' % 'RL GN WM' 
'RAPA_19121977' % 'RL GN WM'  % sp1
'REOO_10031989' % 'RL GN WM'  % sp1
'RFCR_12111978' % 'RL GN WM'  % sp1
'RSVT_17091955' % 'RL GN WM'  % sp1
'SAAA_10051991' % 'RL GN WM'  % sp1
'SAJC_22021972' % 'RL GN WM'  % sp1
'SAVS_23061979' % 'RL GN WM'  % sp1
'UENM_22121979' % 'RL GN WM'  % sp1
'VAPH_11101985' % 'RL GN WM'  % sp1
}' ;


OVERWRITE=true;

PATH_MAT =  '/Volumes/DBNC_03/neuroCOVID/DATA/%S/EEG/%D_pro/';
PATH_EEG =  '/Volumes/DBNC_03/neuroCOVID/DATA/%S/EEG/';


CF=[];
CF.type =  'Morlet';
CF.fwin  = [5] ;% cycles per time window for the 'type' MultiTaper ('MTaper')
CF.step  = [10] ; %numero de puntos entre ventanas e.g. = 10
CF.foi  = [1:0.25:5 5.5:0.5:12 13:35]; %[f1f2 f3 f4 f5 f6 ... ] ; eje de frecuencias  e.g.: 
CF.keeptrials = 'file';%  - 'file' - ('file4chan' only for 'Morlet'  type ) 

load chanlocs_65geodesic_HCGSNv10(64)
chanlocs([23,55,61:64]) = [];

for  nS =Sujetos(25:50) %%
for nT = {'WM', 'RL','GN'}
   disp('Sujetos(25:50)')
    PATH_MAT =  '/Volumes/DBNC_03/neuroCOVID/DATA/%S/EEG/%D_pro/';
    PATH_EEG =  '/Volumes/DBNC_03/neuroCOVID/DATA/%S/EEG/';


    
    
    
SU=     nS{1};  
GR = '';
TAREA = nT{1};       %
DATE = '';


PATH_EEG =strrep( strrep(PATH_EEG,'%G',GR) , '%S' , SU);
if isempty(DATE)
   [paso paso2] = unix(['ls  -d ' PATH_EEG '*_pro' ]);
    DATE = paso2(end-12:end-5);
end
PATH_MAT = strrep( strrep( strrep(PATH_MAT,'%G',GR) , '%S' , SU), '%D', DATE);

if ~exist ([ PATH_MAT 'LAN_' TAREA '_EEG_interp.mat' ], 'file')
  fprintf(['No data  ' SU ' \n'])  
  continue
elseif ~OVERWRITE &&  exist ([ PATH_MAT 'LAN_' TAREA '_EEG_interp_freq.mat' ], 'file')
  fprintf(['freq OK ' SU ' \n'])  
  continue
end


load ([ PATH_MAT 'LAN_' TAREA '_EEG_interp' ])
cd([PATH_MAT ])
if OVERWRITE
unix([' rm ' PATH_MAT '*' TAREA '_cond_1_G_*.ldt' ])
end

% re-reference to mean and delete the no-EEG channels
for t = 1:LAN.trials
    if size(LAN.data{t},1)>58
        LAN.data{t}(59:end,:,:) = [];
    end
    LAN.data{t}(59,:,:) = 0;
    LAN.data{t} = LAN.data{t} - repmat( mean( LAN.data{t} ,1)  ,[59,1,1]);
end
    if size(LAN.nbchan)>58
        LAN.tag.mat(59:end,:,:) = [];
    end
LAN.tag.mat(59,:) = 0;
LAN.chanlocs=chanlocs;
LAN = lan_check(LAN);
LAN = freq_lan(LAN,CF);


% save script in the LAN structure
if isempty(mfilename('fullpath'))
name = matlab.desktop.editor.getActiveFilename;
else
name=[mfilename('fullpath') '.m' ];
end
if isfield(LAN,'m_files')
    LAN.m_files(end+1).name  = name;
    LAN.m_files(end).content  = fileread(name);
    LAN.m_files(end).date = date; 
else
    LAN.m_files(1).name  = name;
    LAN.m_files(1).content  = fileread(name);
    LAN.m_files(1).date = date; 
end



save ([ PATH_MAT 'LAN_' TAREA '_EEG_interp_freq' ],'LAN', '-v7.3')
%save(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG/' DATE '_pro/LAN_' TAREA '_EEG_interp_freq' ] , 'LAN', '-v7.3')
clear LAN
end
end

%%
% %% for checking data de data is OK
if 0 
FT = lan_getdatafile(LAN.freq.powspctrm.filename, LAN.freq.powspctrm.path,LAN.freq.powspctrm.trials);
FT = normal_z(mean( cat(4,FT{:}) ,4));
figure
pcolor2(LAN.freq.time, LAN.freq.freq, squeeze(FT(:,25,:))), shading flat
end

%%







SUs ={
'ASCS_31121976'  %;'08012022';
'CACB_27091965'  %;'17122021';
 'BCUV_02091982' %;'10122021';
'C_DG_15111995'  %; '30062021';
'CACB_27091965'
'CAMO_08111988'  %;'29092021';
 'CARM_14011974'  %;  %19032021_pro
 'CDPR_25031986' %;'29072021';
 'CEBT_30091989' %'03112021';
  'CMCO_17051975' %;'06102021';
 'EERV_01041955' %'25082021'
 'ELFN_09021961' %;'09092021';
 'FABN_19111984'
 'FJAM_10012001'
 'GADM_08111983'
 'GAQJ_21061987'
 'IABF_16091997'
 'IABM_03061982'
 'JAPV_03111995'
 'JAZV_27081991'
 'JDRP_06081956'
 'JFFA_04101976'
 'JFSO_05061967'
 'JGRF_24091991'
 'LAAB_29101982'
 'LDEV_31011988'
 'LMRG_21041962'
 'MAEM_17051957'
 'MFCM_21011984'
 'MGRN_22101979'
 'MLAD_09021957'
 'MLCB_11041981'
 'MPMM_22111977'
 'MTVV_27111974'
 'OAPC_01041962'
 'RABG_09091981'
 'RAPA_19121977'
 'REOO_10031989'
 'RFCR_12111978'
 'RSVT_17091955'
 'SAAA_10051991'
 'SAJC_22021972'
 'SAVS_23061979'
 %'UENM_22121979'
 'VAPH_11101985'
};


CF=[];
CF.type =  'Morlet';
CF.fwin  = [5] ;% cycles per time window for the 'type' MultiTaper ('MTaper')
CF.step  = [10] ; %numero de puntos entre ventanas e.g. = 10
CF.foi  = [1:0.25:5 5.5:0.5:12 13:35]; %[f1f2 f3 f4 f5 f6 ... ] ; eje de frecuencias  e.g.: 
CF.keeptrials = 'yes';%  - 'file' - ('file4chan' only for 'Morlet'  type ) 

load chanlocs_65geodesic_HCGSNv10(64)
chanlocs([23,55,61:64]) = [];



for su =1 %1:numel(SUs)
SU = SUs{su}  ;%;

DA = '';

LUGAR = '/Volumes/DBNC_03/neuroCOVID/DATA/';


     if isempty(DA)
        [paso paso2] = unix(['ls  -d '  LUGAR SU '/EEG/*_pro' ]);
         DA = paso2(end-12:end-5);
     end

     
     % re-reference to mean and delete the no-EEG channels
            
            load([ LUGAR SU '/EEG/'  DA '_pro/LAN_GN_EEG_EYE_Baseline.mat'])
     
            LAN = electrode_lan(LAN,'EYE');
            for t = 1:LAN.trials
                if size(LAN.data{t},1)>58
                    LAN.data{t}(59:end,:,:) = [];
                end
                LAN.data{t}(59,:,:) = 0;
                LAN.data{t} = LAN.data{t} - repmat( mean( LAN.data{t} ,1)  ,[59,1,1]);
            end
            
                
            LAN.tag.mat(59,:) = 0;
            LAN.chanlocs=chanlocs;
            LAN = lan_check(LAN);
            LAN = freq_lan(LAN,CF);


            % save script in the LAN structure
            if isempty(mfilename('fullpath'))
            name = matlab.desktop.editor.getActiveFilename;
            else
            name=[mfilename('fullpath') '.m' ];
            end

            save([ LUGAR SU '/EEG/'  DA '_pro/LAN_GN_EEG_EYE_Baseline_freqEEG.mat'],'LAN')
     
     
     
 
     
     
end
