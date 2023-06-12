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
GR = 'pacientes';


for ns =  {
        'KSDG_20051983'
}'% 'RL GN WM'  % sp1; 
    %'EERV_01041955','GAQJ_21061987','JAZV_27081991','JDRP_06081956','JGRF_24091991',
    %'MLCB_11041981'   'MPMM_22111977','OAPC_01041962','REOO_10031989','UENM_22121979'
    %'LAAB_24101982', 'LDEV_31011988','MLAD_09021957' 'JCZD_03071963', 'JFFA_04101976',
    %'JFSO_05061967''FABN_19111984','IABF_16091997','JAPV_03111995' 'CARM_14011974',%
    %'CDPR_25031986','COGC_18022000','ELFN_09021961','IABM_03061982','CAMO_08111988'
SU=ns{1};
DATA = [];

switch SU
    case 'KSDG_20051983'
        no_ica = [9 12 38 58]; %23 55 61:64           
    otherwise 
        no_ica = []; %23 55 61:64
end


for nm = {'GN'};%'RL','RL','WM'
TASK = nm{1};%'RL'


% directories %S = subject; %G=groups (tie of evaluation)
%PATH_LOG =  '/Volumes/DB_neuroCICS_02/datos_investigadores/neuroCOVID/%G/%S/LOG/';
%PATH_TOBI =  '/Volumes/DB_neuroCICS_02/datos_investigadores/neuroCOVID/%G/%S/EYE/';
%PATH_MAT =  '/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/%S/EEG/%D_pro/';
PATH_MAT =  '/Volumes/Alehermosa/TESIS/%G/%S/EEG/%D_pro/';
%PATH_EEG =  '/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/%S/EEG/';
PATH_EEG =  '/Volumes/Alehermosa/TESIS/%G/%S/EEG/';
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
%continue 
end


load ([ PATH_MAT 'LAN_'  TASK '_array_EYE.mat' ])

% automatic detencion of artifact 

% fix emty trials 

for t = find(isemptycell(LAN.data))
    LAN.data{t} = zeros(size(LAN.data{find(~isemptycell(LAN.data),1)}));  
    disp(['>>>>>>  empty trial '  num2str(t) '!!!!'])
end

LAN = lan_check(LAN);
%

% electrode del EEG 
elec = 1:64;

LAN = vol_thr_lan(LAN,150,0,'bad:V',elec);

% LAN = vol_thr_lan(LAN,3,'z','bad:V',elec); % no ocupar aun esta opcion ,
%                                            % hacer su propio scrits

    %           z-score Porcetaje del scpectro 
    cfga.thr       =   [2 0.25] ;%       %   (sd %spectro)
    cfga.tagname   =   'bad:A';%
    cfga.frange    =   [1 100];%
    %cfga.cat =1;%
    cfga.method    =   'f';%'f';% orLAN
    cfga.nch       =   elec;%
    
LAN = fftamp_thr_lan(LAN,cfga);



ica_sel = 1:64;
icasel = ica_sel;
icasel(no_ica) = [];


% paso opcional 
% 
 Frontales = [1 5 29 32 33 34 37 63 64];
 
 exclude_ele = unique([Frontales no_ica]);

 ica_sel(exclude_ele) = [];% frontales 

    %%% aseguarra sufienctes trials para el ICA
    n=1; 
    LAN.accept = sum(LAN.tag.mat(ica_sel,:),1)<=n;
    while sum(LAN.accept)<(LAN.trials*0.5);
        n=n+1;
        LAN.accept = sum(LAN.tag.mat(ica_sel,:),1)<=n;
    end
    
    
    disp(' ')
    disp('******************************')
    disp([' Ensayos para el ICA::  ' num2str(sum( LAN.accept)) '   '])
    disp('******************************')
    
    
    data1 = cat(2,LAN.data{logical(LAN.accept)});
    data2 = data1(icasel,:);

     % se puede filtrar el dato para el ica en sun pasa alto 1Hz para
     % mejor la deteccion d elso pestañeos 

     [weights,sphere] = runica(     data2    ,'extended', 1,'pca', min(fix(numel(icasel)*0.95),numel(icasel)-1))  ;
     LAN.ica_weights = weights;
     LAN.ica_sphere = sphere;
     LAN.ica_select = icasel;

     
     
load chanlocs_64_Rnet_brainamp.mat
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