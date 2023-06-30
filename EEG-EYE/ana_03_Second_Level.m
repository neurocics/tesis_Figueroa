clear
Sujetos={{
    'CFOI_23101987' % 'RL WM GN'     %
'CAMS_24101973' % 'RL WM GN'  %
'FLFS_08081990' % 'RL GN WM'  % 
'GAMY_03071978' % RL GN WM'  % 
'K_BA_22111991' % RL GN WM % 
'KSDG_20051983' % RL GN WM %  
'MMPM_06011977' % 'RL GN WM' 
'S_WL_25051979' % 'RL GN WM'
'CAPG_27061969' % 'RL GN WM'
%'MJDD_01101984'
}',{% 'RL GN WM' 
'CEAS_23071992' % 'RL GN WM' 
'CAHG_27061988' % 'RL GN WM' 
'YDCL_11021994' % 'RL GN WM' 
'B_FE_24011989' % 'RL GN WM' 
'MEFR_14061991' % 'RL GN WM' 
%'SNCS_22121989' % 'RL GN WM' 
'MAAO_15081988' % 'RL GN WM'  
'GRCN_21111979' % 'RL GN WM'  
'KXMS_13061971' % 'RL GN WM' 
'CFPL_08021984' % 'RL GN WM' 
}'} ;

GRS={'pacientes','controles'};
 

for nG = 1:2
    GR= GRS{nG};
    Sujetos_g = Sujetos{nG};
    numS=0;
for  nS =Sujetos_g %%
       numS=numS+1;
       SU=     nS{1};  
       PATH_EEG =  '/Volumes/Alehermosa/TESIS/%G/%S/EEG/';
       PATH_EEG =strrep( strrep(PATH_EEG,'%G',GR) , '%S' , SU);
       [paso paso2] = unix(['ls  -d ' PATH_EEG '*_pro' ]);
       DATE = paso2(end-12:end-5);      
       Sujetos_New{nG}{numS} = [ SU '/EEG/' DATE '_pro' ];

end
end








%% WM

cfg=[];
cfg.subject  = Sujetos_New;
cfg.groupname =  GRS;

% %   conditionname = {'task', 'rest', ...}
cfg.    comp 	=[2 3 4 5 6  2 3 4 5 6  ] 
                 %[2 3  2 3 ];  	% INDEX OF THE CONDITION TO COMPARED
cfg.    group   = [1 1 1 1 1  2 2 2 2 2  ];
                  %[1 1  2 2 ];  % INDEX OF THE CONDITION TO COMPARED
%     to compare two contrast
%       comp = [ c1 c2 c3 c4]
%       group = [g1 g1 g2 g2]
%       matdif   = [1 -1 2 -2]
%       matdif_transform = 'log', 'none', 'log10'
%       
 cfg.stata=0 ;  
 cfg.paired=0;
% cfg.glm_matrix = [ones(size(Sujetos))'];
% cfg.glm_matrix = [ANOSMIA HOSPITAL Tiempo];
% cfg.RegressorOI=1;
 
 cfg.data_type = 't';       
%    alpha 	=0.05;
%    m		='d'; OR ='i' 	% RELATIONSHEAP TO THE SAMPLES 'i'NDEPENDENT OR 'd'EPENDET
%    bl		=[ 0 0.4];	% BASELINE
%    norma   = 'mdB'
cfg.    mcp      = 0; %or 0
cfg.    nrandom  = 1000;
%    delelectrode = [elec_1 elec_n ... ] % ELECTRODES EXCLUIDED TO THE ANALISIS
%cfg.    mean_eoi = [59] %Electrodos of interes for mean analsis 
%    savesub = 0  % save matrix per subjects  
%cfg.    range = [1 10];
%cfg.time = [-.5 1];
%    permute_name = { 'Real', 'Permut1', 'Permutex', ...}
cfg.matname  = 'LAN';       %                          special carater:
%   .................                                     %S subjectname
cfg.filename  = '/Volumes/Alehermosa/TESIS/%G/%S/LAN_WM_MODEL_ML_SMP_MLxSMP_noise2';
                %'/Volumes/Alehermosa/TESIS/%G/%S/LAN_GN_MODEL_NGO_nGO_nNGP_G_pG_Lf_noise2';
                %'/Volumes/Alehermosa/TESIS/%G/%S/LAN_WM_MODEL_ML_SMP';
                    %G groupname 
                                                         %C conditionname  
                                                        %R permute indicator
%cfg.mcp_fast = 0.5 ;                                                      
GLAN = timefreq_stata([],cfg);

%save 'MS_WM_MODEL_ML_SMP_noise2';
save WM_M4_P10C10

%% GN

cfg=[];
cfg.subject  = Sujetos_New;
cfg.groupname =  GRS;

% %   conditionname = {'task', 'rest', ...}
cfg.    comp 	=[2 3 4 5 6  2 3 4 5 6  ] 
                 %[2 3  2 3 ];  	% INDEX OF THE CONDITION TO COMPARED
cfg.    group   = [1 1 1 1 1  2 2 2 2 2  ];
                  %[1 1  2 2 ];  % INDEX OF THE CONDITION TO COMPARED
%     to compare two contrast
%       comp = [ c1 c2 c3 c4]
%       group = [g1 g1 g2 g2]
%       matdif   = [1 -1 2 -2]
%       matdif_transform = 'log', 'none', 'log10'
%       
 cfg.stata=0 ;  
 cfg.paired=0;
% cfg.glm_matrix = [ones(size(Sujetos))'];
% cfg.glm_matrix = [ANOSMIA HOSPITAL Tiempo];
% cfg.RegressorOI=1;
 
 cfg.data_type = 't';       
%    alpha 	=0.05;
%    m		='d'; OR ='i' 	% RELATIONSHEAP TO THE SAMPLES 'i'NDEPENDENT OR 'd'EPENDET
%    bl		=[ 0 0.4];	% BASELINE
%    norma   = 'mdB'
cfg.    mcp      = 0; %or 0
cfg.    nrandom  = 1000;
%    delelectrode = [elec_1 elec_n ... ] % ELECTRODES EXCLUIDED TO THE ANALISIS
%cfg.    mean_eoi = [59] %Electrodos of interes for mean analsis 
%    savesub = 0  % save matrix per subjects  
%cfg.    range = [1 10];
%cfg.time = [-.5 1];
%    permute_name = { 'Real', 'Permut1', 'Permutex', ...}
cfg.matname  = 'LAN';       %                          special carater:
%   .................                                     %S subjectname
cfg.filename  = '/Volumes/Alehermosa/TESIS/%G/%S/LAN_GN_MODEL_NGO_nGO_nNGP_G_pG_Lf_noise2';
                %'/Volumes/Alehermosa/TESIS/%G/%S/LAN_WM_MODEL_ML_SMP';
                    %G groupname 
                                                         %C conditionname  
                                                        %R permute indicator
%cfg.mcp_fast = 0.5 ;                                                      
GLAN = timefreq_stata([],cfg);

%save 'MS_WM_MODEL_ML_SMP_noise2';
%save WM_M4_P10C10

if 0 % fix time 
load('LAN_WM_MODEL_ML_SMP_noise2.mat')
GLAN.timefreq.time=LAN.freq.model.time;
    
end