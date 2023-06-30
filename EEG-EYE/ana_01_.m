% Time-frequiency calculations 
clear
% eeg analisis 
Sujetos = {{
'CFOI_23101987' % 'RL WM GN'     %
'CAMS_24101973' % 'RL WM GN'  %
'FLFS_08081990' % 'RL GN WM'  % 
'GAMY_03071978' % RL GN WM'  % 
'K_BA_22111991' % RL GN WM % 
'KSDG_20051983' % RL GN WM %  
'MMPM_06011977' % 'RL GN WM' 
'S_WL_25051979' % 'RL GN WM'
'CAPG_27061969' % 'RL GN WM'
'MJDD_01101984'
}',{% 'RL GN WM' 
'CEAS_23071992' % 'RL GN WM' 
'CAHG_27061988' % 'RL GN WM' 
'YDCL_11021994' % 'RL GN WM' 
'B_FE_24011989' % 'RL GN WM' 
'MEFR_14061991' % 'RL GN WM' 
'SNCS_22121989' % 'RL GN WM' 
'MAAO_15081988' % 'RL GN WM'  
'GRCN_21111979' % 'RL GN WM'  
'KXMS_13061971' % 'RL GN WM' 
'CFPL_08021984' % 'RL GN WM' 
}'} ;


GRS={'pacientes','controles'};

OVERWRITE=false;



CF=[];
CF.type =  'Morlet';
CF.fwin  = [5] ;% cycles per time window for the 'type' MultiTaper ('MTaper')
CF.step  = [10] ; %numero de puntos entre ventanas e.g. = 10
CF.foi  = [1:0.25:5 5.5:0.5:12 13:35]; %[f1f2 f3 f4 f5 f6 ... ] ; eje de frecuencias  e.g.: 
CF.keeptrials = 'file';%  - 'file' - ('file4chan' only for 'Morlet'  type ) 



% 
load chanlocs_65_Rnet_brainamp.mat

% chanlocs([23,55,61:64]) = [];
for nG = 1:2
    GR= GRS{nG};
    Sujetos_g = Sujetos{nG};
for  nS =Sujetos_g %%
    
SU=     nS{1};  
DATE = '';    
PATH_MAT =  '/Volumes/Alehermosa/TESIS/%G/%S/EEG/%D_pro/';
PATH_EEG =  '/Volumes/Alehermosa/TESIS/%G/%S/EEG/';  
PATH_EEG =strrep( strrep(PATH_EEG,'%G',GR) , '%S' , SU);

if isempty(DATE)
   [paso paso2] = unix(['ls  -d ' PATH_EEG '*_pro' ]);
    DATE = paso2(end-12:end-5);
end
PATH_MAT = strrep( strrep( strrep(PATH_MAT,'%G',GR) , '%S' , SU), '%D', DATE);   
    


for nT = {'GN', 'WM', 'RL'}
 TAREA = nT{1};   
    

 disp(['Sujeto  '  SU ' tarea '  TAREA ])
 
 switch TAREA
    case 'GN'
        CF.toi = -0.5:0.01:1; %disp('-1,5:1.5')
    case 'WM'
        CF.toi = -0.5:0.01:4.8; %disp('-1:5')
    case 'RL'
        CF.toi = -1:0.01:3.5; %disp('-1:4')
end

    
 
 
 



if ~exist ([ PATH_MAT 'LAN_' TAREA '_EEG_interp.mat' ], 'file')
  fprintf(['No data  ' SU ' \n'])  
  continue
elseif ~OVERWRITE &&  exist ([ PATH_MAT 'LAN_' TAREA '_EEG_interp_freq_lapla.mat' ], 'file')
  fprintf(['freq OK ' SU ' \n'])  
  continue
end


load ([ PATH_MAT 'LAN_' TAREA '_EEG_interp' ])
cd([PATH_MAT ])
if OVERWRITE
unix([' rm ' PATH_MAT '*' TAREA '*_cond*.ldt' ])
end

% re-reference to mean and delete the no-EEG channels
for t = 1:LAN.trials
    if size(LAN.data{t},1)>64
        LAN.data{t}(65:end,:,:) = []; % delete EYE TRACKER channels
        LAN.tag.mat(65:end,:) =[] ;
    end
    LAN.data{t}(65,:,:) = 0; % ref
    LAN.data{t} = LAN.data{t} - repmat( mean( LAN.data{t} ,1)  ,[65,1,1]);
    LAN.tag.mat(65,:) = 0;
end

LAN.chanlocs=chanlocs; 
LAN = lan_check(LAN);
Lcfg=[];
% v.0.0.4
% Compute laplace tranformation using CSD toolbox, 
%         see http://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox/
if exist('H','var')
Lcfg.H         = H;%precalulate H matrix
Lcfg.G         = G;%precalulate G matrix
end
% cfg.head      = heade radius =10
% cfg.lambda    = smooth = 0.00001
Lcfg.centred   = 'Cz' ;   
[ LAN H G ]= lan_laplace(LAN,Lcfg);  
    
    




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



save ([ PATH_MAT 'LAN_' TAREA '_EEG_interp_freq_lapla' ],'LAN', '-v7.3')
%save(['/Volumes/Alehermosa/TESIS/' SU '/EEG/' DATE '_pro/LAN_' TAREA '_EEG_interp_freq' ] , 'LAN', '-v7.3')
clear LAN
end
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



