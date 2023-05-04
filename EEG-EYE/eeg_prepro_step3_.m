% This is not a fullt automatic scritp
%   several part depent of  visual-humnan selection od artefact!!!! 
% This script DO:  
%               1. select de ocular component 
%                  
% 
% dependency:
% LAN toolbox   <*LAN)<]  >= v.1.9.8
% 
%  7 Septimebre  2021
%  Pablo Billeke 

%% Part 1:  Select de ocular component MANUALY
clear
clc

%  'MPMM_22111977','OAPC_01041962','REOO_10031989','UENM_22121979' 
%  'LAAB_24101982', 'LDEV_31011988','MLAD_09021957' 'JCZD_03071963', 
%  'JFFA_04101976','JFSO_05061967''FABN_19111984','IABF_16091997','JAPV_03111995'
%  'CARM_14011974','CDPR_25031986','COGC_18022000'

SU=        'UENM_22121979';%               'SAJC_22021972';%               % 
GR = '';
TAREA = 'GN';        %{'WM', 'RL','GN'}
DATE = '';


PATH_MAT =  '/Volumes/DBNC_03/neuroCOVID/DATA/%S/EEG/%D_pro/';
%PATH_EEG =  '/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/%S/EEG/';
PATH_EEG =  '/Volumes/DBNC_03/neuroCOVID/DATA/%S/EEG/';

PATH_EEG =strrep( strrep(PATH_EEG,'%G',GR) , '%S' , SU);
if isempty(DATE)
   [paso paso2] = unix(['ls  -d ' PATH_EEG '*_pro' ]);
    DATE = paso2(end-12:end-5);
end


PATH_MAT = strrep( strrep( strrep(PATH_MAT,'%G',GR) , '%S' , SU), '%D', DATE);

if exist ([ PATH_MAT 'LAN_' TAREA '_EEG_interp' ], 'file')
  fprintf(['\n Already OK  \n\n'])  
else
load ([ PATH_MAT 'LAN_'  TAREA '_array_EYE_ica.mat' ])

fprintf(['\nFind out ocular and heart-beat components \n\n'])
prepro_plot(LAN) % select and delect ICAS
end

%%
fprintf(['\ndeleted component index: ' num2str(LAN.ica_del) '\n\n'])

LAN.accept(:)=true;
try
LAN = rmfield(LAN,'tag');
end

if strcmp(TAREA,'RL')
   if LAN.time(1,2) >  4.0010
       
        cfg             = [];
        cfg. source     = 'RT'; 
        cfg. ref        = unique(LAN.RT.est);  %[S010 S015 S030 S031] ;   % code of event mark of epoch (references)
        cfg. epoch       = true;
        cfg. times       = [-1.5 4]; % % time for segmentation 
        LAN = lan_latency(LAN, cfg );
       
   end
end


% psoiction chnage in some case 

switch SU  %%                     LABEL --> POSITION ;
    case 'AGBB_26121972'
         elec_shift = [];%s.c.
    case 'BBMA_19041977'
        elec_shift = [63  20;  64 22  ; 62  37; 55 50];
        DATA = '19112021';  
    case 'BCUV_02091982'
        elec_shift = [63  20;  64 22  ; 62  37; 55 50];
        DATA = '10122021';          
    case 'C_DG_15111995'
         elec_shift = [];%s.c.
    case 'CACB_27091965'
         elec_shift = [];%s.c.  
        DATA = '17122021';          
    case 'CARM_14011974'
         elec_shift = []; 
    case 'CDPR_25031986' % revisar, parese que solo para WM !!!!!!!!  63 61 ???
         elec_shift = [63  32;  61   43  ;  23  14  ; 62  60; 55 50];
         if strcmp(TAREA,'WM')
            elec_shift = [63  32;  61   43  ;  23  14  ; 62  60; 55 50];
         end
    case 'CEBT_30091989'
        elec_shift = [55 56; 61 58; 63 16; 64 6; 23 15];
        DATA = '03112021'; 
    case 'CJGV_24021957'
        elec_shift = [64 22; 63 20; 62 37];% sin cambios
        DATA = '20102021';           
    case  'COGC_18022000'
         elec_shift = [64   22; 55   50;  63   20];    
    case 'ELFN_09021961';
          elec_shift = [64 22; 55 50; 63 20; 62 37];     
    case  'FABN_19111984'
         elec_shift = [23   14; 63 32; 55 50; 62  60];
         if strcmp(TAREA,'WM')
             elec_shift = [elec_shift; 64  28];
         end     
    case 'FJAM_10012001'
        elec_shift = [23 12; 63 13];% sin cambios
        DATA = '26112021';              
     case 'FRCA_12081956'
        elec_shift = [];% sin cambios
        DATA = '07012022';            
     case 'GADM_08111983'
        elec_shift = [64 6; 23 15; 61 58; 55 56; 63 16];% sin cambios
        DATA = '26112021';        
    case 'GAQJ_21061987'
         elec_shift = [64 22; 55 50];
    case  'IABF_16091997'
         elec_shift = [];
    case 'IABM_03061982'
         elec_shift = [];% s.c., 
         DATA = '19082021';             
    case 'JAZV_27082021'
         elec_shift = [];%s.c., 
    case  'JAPV_03111995'
         elec_shift = [23    15 ;  55   56 ; 60    58];  
    case  'JCZD_03071963'
         elec_shift = []; %s.c.  
    case 'JDRP_06081956'
         elec_shift =[ 63  32;  61  43 ; 23  14; 62  60; 55 50];
    case  'JFFA_04101976'
         elec_shift = []; %s.c.
    case  'JFSO_05061967'
         elec_shift = [];%s.c.
    case  'JGRF_24091991'
         elec_shift = [];%s.c.
    case   'LAAB_24101982'
         elec_shift = [ ];%s.c.
    case   'LDEV_31011988'
         elec_shift = [];  % s.c.   
     case   'LMRG_21041962'
         elec_shift = [];  % s.c.  
     case 'MAEM_17051957'
        elec_shift = [];  % s.c.  
        DATA = '15092021';  
    case 'MFCM_21011984'
        elec_shift = [23 15; 55 56; 61  58; 63 9];
        DATA = '20102021';   
     case 'MGRN_22101979'
        elec_shift = [];  % s.c. 
        DATA = '08012022';                  
    case 'MTVV_27111974'
        elec_shift = [];  % s.c. 
        DATA = '15092021';     
    case   'MLAD_09021957'
         elec_shift = [23   14;  55  50 ; 62   60];
    case 'MLCB_11041981'
         elec_shift = [64 22;  55  50;  62   37;  63  20];
    case   'MPMM_22111977'
         elec_shift = [];   
    case   'OAPC_01041962'
         elec_shift = []; %s.c, 
     case 'RAPA_19121977'
         elec_shift =  [64 6; 23 15; 61 58; 55 56; 63 16];
        DATA = '03122021';    
     case 'RFCR_12111978'
         elec_shift = []; %s.c, 
        DATA = '07012022';  
    case   'REOO_10031989'
         elec_shift = []; % s.c.        
     case 'RSVT_17091955'
        elec_shift = []; % s.c.
        DATA = '22092021';     
     case 'SAAA_10051991'
        elec_shift = [55 56; 61 58; 64 6; 23 15; 63 16];
        DATA = '07012022';             
     case 'SAJC_22021972'
         elec_shift = []; % s.c.
        DATA = '22092021';  
    case 'SAVS_2306197'
        elec_shift =  [63 16; 23 15; 61 58; 55 56; 64 6];
        DATA = '27102021';                    
    case   'UENM_22121979'
         elec_shift = [ ]; %s.c. 
    case   'CMCO_17051975' 
         elec_shift = [ ]; %s.c. 
    case 'RABG_09091981'     
         %no_ica = [12 13  61 62 64 ]; 
         elec_shift = [63 13;  23  12];
    case 'VAPH_11101985'
        elec_shift =  [64 6; 23 15; 63 16; 55 56; 61 58];
        DATA = '17122021';              
    otherwise 
        elec_shift = [];
end

if ~ isempty(elec_shift)
    fprintf(['\nRemplazos    POCISION <---- LABEL \n'])
for n =  1:size(elec_shift,1)  
    fprintf(['\n                ' num2str(elec_shift(n,2)) '    <----  ' num2str(elec_shift(n,1)) ' \n'])
   for t=1:LAN.trials 
       % REMPLAZOS !!!!
       %                 posicion   <---  label    !!
       LAN.data{t}(elec_shift(n,2),:)=LAN.data{t}(elec_shift(n,1),:);
       % delete repeated channels 
       LAN.data{t}(elec_shift(n,1),:)=0;
   end
end
end  
% elimina electrode en la cara Y OREJAS 
LAN = electrode_lan(LAN,{ 'E64' , 'E63' ,  'E62' ,  'E61' , 'E55' , 'E23'});



%%% REVISAR LOS TRALAS ELIMINADOS Y CANALES MARCADOS 

elec = 1:58;


	if 1
	    d1 = designfilt('lowpassiir','FilterOrder',8, ...
	        'HalfPowerFrequency',47.5,'DesignMethod','butter', 'SampleRate',LAN.srate);
	     %%%% SOLO sirve para el continuo !!!
	    for t=1:LAN.trials    
	    LAN.data{t}(elec,:) = single(filtfilt(d1,double(LAN.data{t}(elec,:)')))';
	    end
    end


    


LAN.accept(:)=1;
LAN = rmfield(LAN, 'tag');
LAN = lan_check(LAN);
LAN = vol_thr_lan(LAN,150,'bad:V',elec);
    cfga.thr    =    [2 0.2] ;%       %   (sd %spectro)
    cfga.tagname=    'bad:A';%
    cfga.frange=     [1 45];%
    %cfga.cat =1;%
    cfga.method =  'f';%'f';% orLAN
    cfga          .nch = elec;%
LAN = fftamp_thr_lan(LAN,cfga);  

id = numel(LAN.tag.labels)+1;
LAN.tag.labels{id} = 'bad:off';


% Mark turn off electrodes 
for t = 1:LAN.trials
   off_e =  nanmean(abs(LAN.data{t}),2)<0.04; 
   off_e(59:end) = false;  
   LAN.tag.mat(off_e,t)=id; 
end



% marcarmalos todos los que tengan nn  canales detectados 
LAN.accept = sum(LAN.tag.mat(:,:)>0,1)<=15;
fprintf(['\nCary out visual impection \n\n'])
prepro_plot(LAN) % select trials / electrode with remaninde artefact for interpolation 

   %%
% Mark turn off electrodes again post visula impections 
for t = 1:LAN.trials
   off_e =  nanmean(abs(LAN.data{t}),2)<0.04; 
   off_e(59:end) = false;  
   LAN.tag.mat(off_e,t)=id; 
end


cfg=[];
cfg.type='bad';
LAN = lan_interp(LAN,cfg);
%prepro_plot(LAN)


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

fprintf(['\nSaving ... \n\n'])
fprintf(['' PATH_MAT 'LAN_' TAREA '_EEG_interp\n'])
% save in local disck
save([ PATH_MAT '/LAN_' TAREA '_EEG_interp'],'LAN', '-v7.3');

fprintf(['${NUBE}/LAN_' TAREA '_EEG_interp\n'])
%save in neuroCICS DB
%mkdir(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG' ] , [  DATA '_pro' ])
save(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG/' DATE '_pro/LAN_' TAREA '_EEG_interp' ] , 'LAN', '-v7.3')

%clear
fprintf(['\nDone \n\n'])

