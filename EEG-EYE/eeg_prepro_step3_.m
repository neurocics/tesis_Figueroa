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



SU=        'KSDG_20051983';%               'SAJC_22021972';%               % 
GR = 'pacientes';
TAREA = 'GN';        %{'WM', 'RL','GN'}
DATE = '';


PATH_MAT =  '/Volumes/Alehermosa/TESIS/%G/%S/EEG/%D_pro/';
%PATH_EEG =  '/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/%S/EEG/';
PATH_EEG =  '/Volumes/Alehermosa/TESIS/%G/%S/EEG/';

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

%% Part 2
fprintf(['\ndeleted component index: ' num2str(LAN.ica_del) '\n\n'])

LAN.accept(:)=true;
try
LAN = rmfield(LAN,'tag');
LAN = lan_check(LAN);
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



 
% elimina electrode en la cara Y OREJAS 
%LAN = electrode_lan(LAN,{ 'E64' , 'E63' ,  'E62' ,  'E61' , 'E55' , 'E23'});



%%% REVISAR LOS TRALAS ELIMINADOS Y CANALES MARCADOS 

elec = 1:64;


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
                                                        
LAN = vol_thr_lan(LAN,150,'v','bad:V',elec);
LAN = vol_thr_lan(LAN,4,'z','bad:Vz',elec);
LAN = vol_thr_lan(LAN,3.5,'c','bad:Vc',elec);

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
   off_e(65:end) = false;  
   LAN.tag.mat(off_e,t)=id; 
end



% marcarmalos todos los que tengan nn  canales detectados 
LAN.accept = sum(LAN.tag.mat(:,:)>0,1)<=15; % se cambió a 10 canales malos por trial 
fprintf(['\nCary out visual impection \n\n'])
prepro_plot(LAN) % select trials / electrode with remaninde artefact for interpolation 

   %% Part 3
% Mark turn off electrodes again post visula impections 
for t = 1:LAN.trials
   off_e =  nanmean(abs(LAN.data{t}),2)<0.04; 
   off_e(65:end) = false;  
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
%save(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG/' DATE '_pro/LAN_' TAREA '_EEG_interp' ] , 'LAN', '-v7.3')

%clear
fprintf(['\nDone \n\n'])

