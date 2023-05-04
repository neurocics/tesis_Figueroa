% clear
% % eeg analisis 
% Sujetos = {
% 'AGBB_26121972' % 'RL WM GN'     % OK
% 'ASCS_31121976' % 'RL WM GN'  % sp1
% 'BBMA_19041977' % 'RL GN WM'  % sp1
% 'BCUV_02091982' % RL GN WM'  % sp1
% 'C_DG_15111995' % RL GN WM % sp1 GN
% 'CACB_27091965' % RL GN WM % sp1 
% 'CAMO_08111988' % 'RL GN WM' 
% 'CARM_14011974' % 'RL GN WM'
% 'CDPR_25031986' % 'RL GN WM'
% 'CEBT_30091989' % 'RL GN WM' % sp1 
% 'CJGV_24021957' % 'RL GN WM' % sp1 
% 'CMCO_17051975' % 'RL GN WM' 
% 'COGC_18022000' % 'RL GN WM' 
% 'EERV_01041955' % 'RL GN WM' 
% 'ELFN_09021961' % 'RL GN WM' 
% 'FABN_19111984' % 'RL GN WM' 
% 'FJAM_10012001' % 'RL GN WM' % sp1 
% 'FRCA_12081956' % 'RL GN WM' % sp1 
% 'GADM_08111983' % 'RL GN WM' % sp1 
% 'GAQJ_21061987' % 'RL GN WM' 
% 'IABF_16091997' % 'RL GN WM' 
% 'IABM_03061982' % 'RL GN WM' %RRR
% 'JAPV_03111995' % 'RL GN WM'
% 'JAZV_27081991' % 'RL GN WM'
% 'JCZD_03071963' % 'RL GN WM'% RRR % sp1 
% 'JDRP_06081956' % 'RL GN WM'  % sp1 
% 'JFFA_04101976' % 'RL GN WM'  % sp1 
% 'JFSO_05061967' % 'RL GN WM'  % sp1 
% 'JGRF_24091991' % 'RL GN WM'  % sp1 
% 'LAAB_29101982' % 'RL GN WM'  % sp1
% 'LDEV_31011988' % 'RL GN WM'  % sp1
% 'LMRG_21041962' % 'RL GN WM'  % sp1
% 'MAEM_17051957' % 'RL GN WM'  % sp1
% 'MFCM_21011984' % 'RL GN WM'  % sp1
% 'MGRN_22101979' % 'RL GN WM'  % sp1
% 'MLAD_09021957' % 'RL GN WM'  % sp1
% 'MLCB_11041981' % 'RL GN WM'  % sp1
% 'MPMM_22111977' % 'RL GN WM'  % sp1
% 'MTVV_27111974' % 'RL GN WM'  % sp1
% 'OAPC_01041962' % 'RL GN WM'  % sp1
% 'RABG_09091981' % 'RL GN WM' 
% 'RAPA_19121977' % 'RL GN WM'  % sp1
% 'REOO_10031989' % 'RL GN WM'  % sp1
% 'RFCR_12111978' % 'RL GN WM'  % sp1
% 'RSVT_17091955' % 'RL GN WM'  % sp1
% 'SAAA_10051991' % 'RL GN WM'  % sp1
% 'SAJC_22021972' % 'RL GN WM'  % sp1
% 'SAVS_23061979' % 'RL GN WM'  % sp1
% 'UENM_22121979' % 'RL GN WM'  % sp1
% 'VAPH_11101985' % 'RL GN WM'  % sp1
% }' ;
% 
% 
%     
%     
%     fprintf(['\n'])
%     for nS=Sujetos;
%     PATH_MAT =  '/Volumes/DBNC_03/neuroCOVID/DATA/%S/EEG/%D_pro/';
%     PATH_EEG =  '/Volumes/DBNC_03/neuroCOVID/DATA/%S/EEG/';
%     PATH_LOG =  '/Volumes/DBNC_03/neuroCOVID/DATA/%S/LOG/';
%     SU=     nS{1};  
%     GR = '';
%     %TAREA = nT{1};       %
%     DATE = '';
% 
% 
%     PATH_EEG =strrep( strrep(PATH_EEG,'%G',GR) , '%S' , SU);
%     if isempty(DATE)
%        [paso paso2] = unix(['ls  -d ' PATH_EEG '*_pro' ]);
%         DATE = paso2(end-12:end-5);
%     end
%     PATH_MAT = strrep( strrep( strrep(PATH_MAT,'%G',GR) , '%S' , SU), '%D', DATE);
%     PATH_LOG = strrep( strrep( strrep(PATH_LOG,'%G',GR) , '%S' , SU), '%D', DATE);
%     fprintf([' '''  PATH_MAT(34:end-1)  '''\n'])
%     end
% 
% 
% 
% %%
% 
Sujetos={
 'AGBB_26121972/EEG/19052021_pro'
 'ASCS_31121976/EEG/08012022_pro'
 'BBMA_19041977/EEG/19112021_pro'
 'BCUV_02091982/EEG/10122021_pro'
 'C_DG_15111995/EEG/30062021_pro'
 'CACB_27091965/EEG/17122021_pro'
 'CAMO_08111988/EEG/29092021_pro'
 'CARM_14011974/EEG/19032021_pro'
 'CDPR_25031986/EEG/29072021_pro'
 'CEBT_30091989/EEG/03112021_pro'
 'CJGV_24021957/EEG/20102021_pro'
 'CMCO_17051975/EEG/06102021_pro'
 'COGC_18022000/EEG/19082021_pro'
 'EERV_01041955/EEG/25082021_pro'
 'ELFN_09021961/EEG/09092021_pro'
 'FABN_19111984/EEG/19082021_pro'
 'FJAM_10012001/EEG/26112021_pro'
 'FRCA_12081956/EEG/07012022_pro'
 'GADM_08111983/EEG/26112021_pro'
 'GAQJ_21061987/EEG/16062021_pro'
 'IABF_16091997/EEG/12032021_pro'
 'IABM_03061982/EEG/19082021_pro'
 'JAPV_03111995/EEG/23062021_pro'
 'JAZV_27081991/EEG/05082021_pro'
 'JCZD_03071963/EEG/24022021_pro'
 'JDRP_06081956/EEG/05082021_pro'
 'JFFA_04101976/EEG/15072021_pro'
 'JFSO_05061967/EEG/29072021_pro'
 'JGRF_24091991/EEG/24032021_pro'
 'LAAB_29101982/EEG/02062021_pro'
 'LDEV_31011988/EEG/09062021_pro'
 'LMRG_21041962/EEG/08092021_pro'
 'MAEM_17051957/EEG/15092021_pro'
 'MFCM_21011984/EEG/20102021_pro'
 'MGRN_22101979/EEG/08012022_pro'
 'MLAD_09021957/EEG/07072021_pro'
 'MLCB_11041981/EEG/25082021_pro'
 'MPMM_22111977/EEG/19032021_pro'
 'MTVV_27111974/EEG/15092021_pro'
 'OAPC_01041962/EEG/23062021_pro'
 'RABG_09091981/EEG/06102021_pro'
 'RAPA_19121977/EEG/03122021_pro'
 'REOO_10031989/EEG/23062021_pro'
 'RFCR_12111978/EEG/07012022_pro'
 'RSVT_17091955/EEG/22092021_pro'
 'SAAA_10051991/EEG/07012022_pro'
 'SAJC_22021972/EEG/22092021_pro'
 'SAVS_23061979/EEG/27102021_pro'
 'UENM_22121979/EEG/19052021_pro'
 'VAPH_11101985/EEG/17122021_pro'
    }';


ANOSMIA=[
0
1
1
0
1
0
1
0
0
0
0
0
1
1
1
1
1
1
1
1
1
1
1
1
1
1
1
0
1
1
0
1
1
1
1
0
0
1
0
1
1
0
1
0
0
0
1
1
0
0];


HOSPITAL=[
0
1
0
0
0
0
1
1
0
1
1
1
0
0
1
0
0
0
1
0
0
1
0
0
0
0
1
1
0
0
0
1
1
1
0
1
1
1
1
1
1
0
0
1
1
0
1
1
0
0    
];


EDAD=[
49
45
44
39
36
56
33
47
35
32
64
46
21
66
60
37
19
65
38
34
24
39
26
30
57
65
45
55
30
39
33
59
64
37
42
64
40
44
47
59
40
44
32
43
66
30
51
42
42
36
];



Tiempo=[%
    2
5
2
1
11
5
11
8
14
8
13
7
15
6
14
16
3
5
7
4
4
13
12
13
7
6
8
13
10
4
12
14
14
6
3
13
11
1
13
8
14
5
6
6
14
9
14
9
2
1   ];




%%
cfg=[];
cfg.subject{1}  =  Sujetos(ANOSMIA==1);
cfg.subject{2}  =  Sujetos(ANOSMIA==0);
cfg.groupname =  {'DATA','DATA'};

% %   conditionname = {'task', 'rest', ...}
cfg.    comp 	=[2 2 ];  	% INDEX OF THE CONDITION TO COMPARED
cfg.    group    =[1 2];  % INDEX OF THE CONDITION TO COMPARED
%     to compare two contrast
%       comp = [ c1 c2 c3 c4]
%       group = [g1 g1 g2 g2]
%       matdif   = [1 -1 2 -2]
%       matdif_transform = 'log', 'none', 'log10'
%       
 cfg.stata=1 ;  
 cfg.paired=0;
% cfg.glm_matrix = [ones(size(Sujetos))'];
% cfg.glm_matrix = [ANOSMIA HOSPITAL Tiempo];
% cfg.RegressorOI=1;
 
 cfg.data_type = 't';       
%    alpha 	=0.05;
%    m		='d'; OR ='i' 	% RELATIONSHEAP TO THE SAMPLES 'i'NDEPENDENT OR 'd'EPENDET
%    bl		=[ 0 0.4];	% BASELINE
%    norma   = 'mdB'
cfg.    mcp      = 1; %or 0
cfg.    nrandom  = 1000;
%    delelectrode = [elec_1 elec_n ... ] % ELECTRODES EXCLUIDED TO THE ANALISIS
%cfg.    mean_eoi = [59] %Electrodos of interes for mean analsis 
%    savesub = 0  % save matrix per subjects  
%cfg.    range = [1 10];
%cfg.time = [-.5 1];
%    permute_name = { 'Real', 'Permut1', 'Permutex', ...}
cfg.matname  = 'LAN';       %                          special carater:
%   .................                                     %S subjectname
cfg.filename  = '/Volumes/DBNC_03/neuroCOVID/DATA/%S/LAN_RL_MODEL_DM_r_Uc_LD_Ri';                                     %G groupname 
                                                         %C conditionname  
                                                        %R permute indicator
cfg.mcp_fast = 0.5 ;                                                      
GLAN = timefreq_stata([],cfg);




save model_AnosmiaWilcoxon_MCC_1000

%%
%%
cfg=[];
cfg.subject{1}  =  Sujetos; %(ANOSMIA==1);
cfg.subject{2}  =  Sujetos(ANOSMIA==0);
cfg.subject{3}  =  Sujetos(ANOSMIA==1);
cfg.subject{4}  =  Sujetos(HOSPITAL==0);
cfg.subject{5}  =  Sujetos(HOSPITAL==1);

cfg.groupname =  {'DATA','DATA','DATA','DATA','DATA'};

% %   conditionname = {'task', 'rest', ...}
cfg.    comp 	=[2  2 2 2 2 ];  	% INDEX OF THE CONDITION TO COMPARED
cfg.    group    =[1 2 3 4 5 ];  % INDEX OF THE CONDITION TO COMPARED
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
cfg.filename  = '/Volumes/DBNC_03/neuroCOVID/DATA/%S/LAN_RL_MODEL_DM_r_Uc_LD_Ri_noise2';                                     %G groupname 
                                                         %C conditionname  
                                                        %R permute indicator
cfg.mcp_fast = 0.5 ;                                                      
GLAN = timefreq_stata([],cfg);




save /Volumes/DBNC_03/neuroCOVID/ANA/EEG/model_all_noise
%%

%%
cfg=[];
cfg.subject{1}  =  Sujetos;%(ANOSMIA==1);
%cfg.subject{2}  =  Sujetos(ANOSMIA==0);
cfg.groupname =  {'DATA'};

% %   conditionname = {'task', 'rest', ...}
cfg.    comp 	=[4  ];  	% INDEX OF THE CONDITION TO COMPARED
cfg.    group    =[1 ];  % INDEX OF THE CONDITION TO COMPARED
%     to compare two contrast
%       comp = [ c1 c2 c3 c4]
%       group = [g1 g1 g2 g2]
%       matdif   = [1 -1 2 -2]
%       matdif_transform = 'log', 'none', 'log10'
%       
 cfg.stata='glm' ;  
 cfg.paired=0;
 %cfg.glm_matrix = [ones(size(Sujetos))'];
 cfg.glm_matrix = [ones(size(Sujetos))' ANOSMIA.*(-1) HOSPITAL Tiempo];
 cfg.RegressorOI=2;
 
 cfg.data_type = 't';       
%    alpha 	=0.05;
%    m		='d'; OR ='i' 	% RELATIONSHEAP TO THE SAMPLES 'i'NDEPENDENT OR 'd'EPENDET
%    bl		=[ 0 0.4];	% BASELINE
%    norma   = 'mdB'
cfg.    mcp      = 1; %or 0
cfg.    nrandom  = 1000;
%    delelectrode = [elec_1 elec_n ... ] % ELECTRODES EXCLUIDED TO THE ANALISIS
%cfg.    mean_eoi = [59] %Electrodos of interes for mean analsis 
%    savesub = 0  % save matrix per subjects  
%cfg.    range = [1 10];
%cfg.time = [-.5 1];
%    permute_name = { 'Real', 'Permut1', 'Permutex', ...}
cfg.matname  = 'LAN';       %                          special carater:
%   .................                                     %S subjectname
cfg.filename  = '/Volumes/DBNC_03/neuroCOVID/DATA/%S/LAN_WM_MODEL_ML_SMP_MLxSMP_noise2';                                     %G groupname 
                                                         %C conditionname  
                                                        %R permute indicator
cfg.mcp_fast = 0.5 ;                                                      
GLAN = timefreq_stata([],cfg);




save model_WM_(SMPxML)_Anosmia_glm_MCC_1000


%%
cfg=[];
cfg.subject{1}  =  Sujetos;%(ANOSMIA==1);
%cfg.subject{2}  =  Sujetos(ANOSMIA==0);
cfg.groupname =  {'DATA'};

% %   conditionname = {'task', 'rest', ...}
cfg.    comp 	=[2  ];  	% INDEX OF THE CONDITION TO COMPARED
cfg.    group    =[1 ];  % INDEX OF THE CONDITION TO COMPARED
%     to compare two contrast
%       comp = [ c1 c2 c3 c4]
%       group = [g1 g1 g2 g2]
%       matdif   = [1 -1 2 -2]
%       matdif_transform = 'log', 'none', 'log10'
%       
 cfg.stata='glm' ;  
 cfg.paired=0;
 %cfg.glm_matrix = [ones(size(Sujetos))'];
 cfg.glm_matrix = [ones(size(Sujetos))' ANOSMIA.*(-1) HOSPITAL Tiempo];
 cfg.RegressorOI=2;
 
 cfg.data_type = 't';       
%    alpha 	=0.05;
%    m		='d'; OR ='i' 	% RELATIONSHEAP TO THE SAMPLES 'i'NDEPENDENT OR 'd'EPENDET
%    bl		=[ 0 0.4];	% BASELINE
%    norma   = 'mdB'
cfg.    mcp      = 1; %or 0
cfg.    nrandom  = 1000;
%    delelectrode = [elec_1 elec_n ... ] % ELECTRODES EXCLUIDED TO THE ANALISIS
%cfg.    mean_eoi = [59] %Electrodos of interes for mean analsis 
%    savesub = 0  % save matrix per subjects  
%cfg.    range = [1 10];
%cfg.time = [-.5 1];
%    permute_name = { 'Real', 'Permut1', 'Permutex', ...}
cfg.matname  = 'LAN';       %                          special carater:
%   .................                                     %S subjectname
cfg.filename  = '/Volumes/DBNC_03/neuroCOVID/DATA/%S/LAN_RL_MODEL_DM_r_Uc_LD_Ri';                                     %G groupname 
                                                         %C conditionname  
                                                        %R permute indicator
cfg.mcp_fast = 0.2 ;                                                      
GLAN = timefreq_stata([],cfg);




save model_RL_(Utility)_Anosmia_glm_MCC_1000

%%


cfg=[];
cfg.subject{1}  =  Sujetos;%(ANOSMIA==1);
%cfg.subject{2}  =  Sujetos(ANOSMIA==0);
cfg.groupname =  {'DATA'};

% %   conditionname = {'task', 'rest', ...}
cfg.    comp 	=[2  ];  	% INDEX OF THE CONDITION TO COMPARED
cfg.    group    =[1 ];  % INDEX OF THE CONDITION TO COMPARED
%     to compare two contrast
%       comp = [ c1 c2 c3 c4]
%       group = [g1 g1 g2 g2]
%       matdif   = [1 -1 2 -2]
%       matdif_transform = 'log', 'none', 'log10'
%       
 cfg.stata='glm' ;  
 cfg.paired=0;
 %cfg.glm_matrix = [ones(size(Sujetos))'];
 cfg.glm_matrix = [ones(size(Sujetos))' ANOSMIA HOSPITAL Tiempo];
 cfg.RegressorOI=2;
 
 cfg.data_type = 't';       
%    alpha 	=0.05;
%    m		='d'; OR ='i' 	% RELATIONSHEAP TO THE SAMPLES 'i'NDEPENDENT OR 'd'EPENDET
%    bl		=[ 0 0.4];	% BASELINE
%    norma   = 'mdB'
cfg.    mcp      = 1; %or 0
cfg.    nrandom  = 500;
%    delelectrode = [elec_1 elec_n ... ] % ELECTRODES EXCLUIDED TO THE ANALISIS
%cfg.    mean_eoi = [59] %Electrodos of interes for mean analsis 
%    savesub = 0  % save matrix per subjects  
%cfg.    range = [1 10];
%cfg.time = [-.5 1];
%    permute_name = { 'Real', 'Permut1', 'Permutex', ...}
cfg.matname  = 'LAN';       %                          special carater:
%   .................                                     %S subjectname
cfg.filename  = '/Volumes/DBNC_03/neuroCOVID/DATA/%S/LAN_GN_MODEL_NGO_nGO_nNGP_G_pG_Lf_noise2';                                     %G groupname 
                                                         %C conditionname  
                                                        %R permute indicator
cfg.mcp_fast = 0.4 ;                                                      
GLAN = timefreq_stata([],cfg);




save model_GN_(NGO)_Anosmia_glm_MCC_500
