% This script DO:   
%               1. caculate  ICA in the 50%-best trails 
%                  
% 
% dependency:
% LAN toolbox   <*LAN)<]  >= v.1.9.8
% 
% 31 Agosto 2021
%  Pablo Billeke 

clear
GR = '';


for ns =  {
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
}'% 'RL GN WM'  % sp1; 
    %'EERV_01041955','GAQJ_21061987','JAZV_27081991','JDRP_06081956','JGRF_24091991',
    %'MLCB_11041981'   'MPMM_22111977','OAPC_01041962','REOO_10031989','UENM_22121979'
    %'LAAB_24101982', 'LDEV_31011988','MLAD_09021957' 'JCZD_03071963', 'JFFA_04101976',
    %'JFSO_05061967''FABN_19111984','IABF_16091997','JAPV_03111995' 'CARM_14011974',%
    %'CDPR_25031986','COGC_18022000','ELFN_09021961','IABM_03061982','CAMO_08111988'
SU=ns{1};
DATA = [];
switch SU
    case 'ASCS_31121976'
        no_ica = [23 55 61:64];% sin cambios
        DATA = '08012022';
    case 'AGBB_26121972'
         no_ica = [15 55 58];
    case 'BBMA_19041977'
        no_ica = [ 23 50 61 22 20 37 50 ];% sin cambios
        DATA = '19112021'; 
    case 'BCUV_02091982'
        no_ica = [ 23 50 61 22 20 37 50 ];% sin cambios
        DATA = '10122021';         
    case 'C_DG_15111995'
         no_ica = [23 55 62 63];
    case 'CACB_27091965'
        no_ica = [ 23 55 61:64];% sin cambios
        DATA = '17122021';       
    case 'CARM_14011974'
         no_ica = [23 32 43 50 53 54 55 60 63];  
    case 'CEBT_30091989'
        no_ica = [56 58 16 6 15   62 ];% sin cambios
        DATA = '03112021';              
    case 'CJGV_24021957'
        no_ica = [23 22 50 20 37 55 61];% sin cambios
        DATA = '20102021';           
    case 'CDPR_25031986'
         no_ica = [14 32 43 45 50 60];      
    case  'COGC_18022000'
         no_ica = [20 22 23 37 50 61]; 
    case  'EERV_01041955'
         no_ica = [ 55 ];  
    case 'ELFN_09021961';
          no_ica = [20 22 23 37 50 61 ];
    case  'GAQJ_21061987'
         no_ica = [ 20 22 23 37 50 61];        
    case  'FABN_19111984'
         no_ica = [14 32 43 50 60 63];            
    case 'FJAM_10012001'
        no_ica = [ 55 61 62 64];% sin cambios
        DATA = '26112021';         
     case 'FRCA_12081956'
        no_ica = [23 55 61:64];% sin cambios
        DATA = '07012022';        
     case 'GADM_08111983'
        no_ica = [ 6 15 58 56 16 62];% sin cambios
        DATA = '26112021';        
    case  'IABF_16091997'
         no_ica = [15];       
    case  'IABM_03061982';
        no_ica = [55]; 
        DATA = '19082021';             
    case  'JAPV_03111995'
         no_ica = [15 56 58 61:64]; 
     case  'JAZV_27081991'
         no_ica = [55 62 ];                
    case  'JCZD_03071963'
         no_ica = [23 62];   
    case  'JDRP_06081956'
         no_ica = [14 28 32 43 50 60 61 ];            
    case  'JFFA_04101976'
         no_ica = [23 55 61:64];
    case  'JFSO_05061967'
         no_ica = [ 47 52 55 63 23 ];
    case  'JGRF_24091991'
         no_ica = [ 15 55 23 ];
    case   'LAAB_24101982'
         no_ica = [ 3 20 22 23 37 50 55 61 ];
    case   'LDEV_31011988'
         no_ica = [ 15 55 58  ]; 
    case   'LMRG_21041962'
         no_ica = [23 55  61:64 ]; 
    case 'MAEM_17051957'
        no_ica = [23 55 61:64];
        DATA = '15092021';     
    case 'MFCM_21011984'
        no_ica = [6 15 56 16 58 9];
        DATA = '20102021';   
    case 'MGRN_22101979'
        no_ica = [23 55 61:64];
        DATA = '08012022';                  
    case   'MLAD_09021957'
         no_ica = [ 14  32 43 50 60];
    case   'MLCB_11041981'
         no_ica = [3 20 22 23 37 50 61 ];         
    case   'MPMM_22111977'
         no_ica = [ 23 55];  
    case 'MGRN_22101979'
        no_ica = [23 55 61:64];
        DATA = '08012022';          
   case 'MTVV_27111974'
        no_ica = [23 55 61:64];
        DATA = '15092021';        
    case   'OAPC_01041962'
         no_ica = [  63];  
     case 'RAPA_19121977'
        no_ica = [6 15 58 56 16 62];
        DATA = '03122021';          
    case   'REOO_10031989'
         no_ica = [ 20 32 43 50 55 60 63]; 
     case 'RFCR_12111978'
        no_ica = [23 55 61:64];
        DATA = '07012022';        
     case 'RSVT_17091955'
        no_ica = [23 55 61:64];
        DATA = '22092021';           
     case 'SAAA_10051991'
        no_ica = [56 58 6 15 16  62];
        DATA = '07012022';        
     case 'SAJC_22021972'
        no_ica = [23 55 61:64];
        DATA = '22092021'; 
     case 'SAVS_2306197'
        no_ica = [16 15 58 56 6 62];
        DATA = '27102021';            
    case   'UENM_22121979'
         no_ica = [ 20 22 23 37 50 55 61];  
    case   'CAMO_08111988'
         no_ica = [ 61:64 ]; 
    case 'CMCO_17051975'     
         no_ica = [ 61:64 ]; 
    case 'RABG_09091981'     
         no_ica = [12 13  61 62 64 ];   
    case 'VAPH_11101985'
        no_ica = [6 15 16 56 58 62];
        DATA = '17122021';               
    otherwise 
        no_ica = [23 55 61:64];
end


for nm = { 'RL','WM','GN'};%
TASK = nm{1};%'RL'


% directories %S = subject; %G=groups (tie of evaluation)
%PATH_LOG =  '/Volumes/DB_neuroCICS_02/datos_investigadores/neuroCOVID/%G/%S/LOG/';
%PATH_TOBI =  '/Volumes/DB_neuroCICS_02/datos_investigadores/neuroCOVID/%G/%S/EYE/';
%PATH_MAT =  '/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/%S/EEG/%D_pro/';
PATH_MAT =  '/Volumes/DBNC_03/neuroCOVID/DATA/%S/EEG/%D_pro/';
%PATH_EEG =  '/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/%S/EEG/';
PATH_EEG =  '/Volumes/DBNC_03/neuroCOVID/DATA/%S/EEG/';
%find DATA

PATH_EEG =strrep( strrep(PATH_EEG,'%G',GR) , '%S' , SU);
if isempty(DATA)
   [paso paso2] = unix(['ls  -d ' PATH_EEG '*_pro' ]);
    DATA = paso2(end-12:end-5);

end

    disp([ 'data of adquision : ' DATA])
% create real PATH for the files 
%PATH_LOG =strrep( strrep(PATH_LOG,'%G',GR) , '%S' , SU);
%PATH_TOBI =strrep( strrep(PATH_TOBI,'%G',GR) , '%S' , SU);
PATH_MAT = strrep(strrep( strrep(PATH_MAT,'%G',GR) , '%S' , SU),'%D',DATA);

if ~exist([ PATH_MAT 'LAN_'  TASK '_array_EYE.mat' ],'file')
disp(['No preprocesing data ' SU '>>' TASK ])    
continue 
end

if exist([ PATH_MAT 'LAN_'  TASK '_array_EYE_ica.mat' ],'file') || exist([ PATH_MAT 'LAN_'  TASK '_EEG_interp.mat' ],'file')
disp(['data already OK  ' SU '>>' TASK ])    
continue 
end


load ([ PATH_MAT 'LAN_'  TASK '_array_EYE.mat' ])

% automatic detencion of artifact 

% fix emty trials 

for t = find(isemptycell(LAN.data))
    LAN.data{t} = zeros(size(LAN.data{find(~isemptycell(LAN.data),1)}));  
    disp(['>>>>>>  empty trail '  num2str(t) '!!!!'])
end
 LAN = lan_check(LAN);
%


elec = 1:64;

LAN = vol_thr_lan(LAN,150,'bad:V',elec);
    cfga.thr    =    [2 0.25] ;%       %   (sd %spectro)
    cfga.tagname=    'bad:A';%
    cfga.frange=     [1 40];%
    %cfga.cat =1;%
    cfga.method =  'f';%'f';% orLAN
    cfga          .nch = elec;%
    
LAN = fftamp_thr_lan(LAN,cfga);



ica_sel = 1:64;
ica_sel(no_ica) = [];
icasel = ica_sel;
ica_sel(ica_sel<14) = [];

    %%% aseguarra sufienctes trials para el ICA
    n=1; 
    LAN.accept = sum(LAN.tag.mat(ica_sel,:),1)<=n;
    while sum(LAN.accept)<(LAN.trials*0.5);
        n=n+1;
        LAN.accept = sum(LAN.tag.mat(ica_sel,:),1)<=n;
    end
    
    
    disp(' ')
    disp('******************************')
    disp([' Esayos para el ICA::  ' num2str(sum( LAN.accept)) '   '])
    disp('******************************')
    
    
    data1 = cat(2,LAN.data{logical(LAN.accept)});
    data2 = data1(icasel,:);

     [weights,sphere] = runica(     data2    ,'extended', 1,'pca',min(fix(numel(icasel)*0.95),numel(icasel)-1));
     LAN.ica_weights = weights;
     LAN.ica_sphere = sphere;
     LAN.ica_select = icasel;

     
     
load chanlocs_65geodesic_HCGSNv10(64)
paso = LAN.chanlocs;
LAN.chanlocs =  chanlocs(1:64);
for e = 65:LAN.nbchan;
   LAN.chanlocs(e).labels = paso(e).labels;
   LAN.chanlocs(e).type = paso(e).type;
end

%%

%prepro_plot(LAN) % select and delect ICAS

%%


if isempty(mfilename('fullpath'))
name = matlab.desktop.editor.getActiveFilename;
else
name=[mfilename('fullpath') '.m' ];
end

%  lanversion >= 1.9.6
%LAN = lan_add_this_m(LAN,name); 
%

if isfield(LAN,'m_files')
    LAN.m_files(end+1).name  = name;
    LAN.m_files(end).content  = fileread(name);
    LAN.m_files(end).date = date; 
else
    LAN.m_files(1).name  = name;
    LAN.m_files(1).content  = fileread(name);
    LAN.m_files(1).date = date; 
end



save ([ PATH_MAT 'LAN_'  TASK '_array_EYE_ica.mat' ], 'LAN', '-v7.3')
clear LAN 
end
end