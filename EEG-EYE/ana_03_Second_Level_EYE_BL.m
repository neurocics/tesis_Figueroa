
clear
clc
% OK 
Sujetos ={
'ASCS_31121976'
'BCUV_02091982'
'C_DG_15111995'
'CACB_27091965'
'CAMO_08111988'
'CARM_14011974'
'CDPR_25031986'
'CEBT_30091989'
'CMCO_17051975'
'EERV_01041955'
'ELFN_09021961'
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
'VAPH_11101985'
};


HOSPITAL=[
1
0
0
0
1
1
0
1
1
0
1
0
0
1
0
0
1
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
];

ANOSMIA=[
 1
0
1
0
1
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
]


    fprintf(['\n'])
    for nS=Sujetos';
    PATH_MAT =  '/Volumes/DBNC_03/neuroCOVID/DATA/%S/EEG/%D_pro/';
    PATH_EEG =  '/Volumes/DBNC_03/neuroCOVID/DATA/%S/EEG/';
    PATH_LOG =  '/Volumes/DBNC_03/neuroCOVID/DATA/%S/LOG/';
    SU=     nS{1};  
    GR = '';
    %TAREA = nT{1};       %
    DATE = '';


    PATH_EEG =strrep( strrep(PATH_EEG,'%G',GR) , '%S' , SU);
    if isempty(DATE)
       [paso paso2] = unix(['ls  -d ' PATH_EEG '*_pro' ]);
        DATE = paso2(end-12:end-5);
    end
    PATH_MAT = strrep( strrep( strrep(PATH_MAT,'%G',GR) , '%S' , SU), '%D', DATE);
    PATH_LOG = strrep( strrep( strrep(PATH_LOG,'%G',GR) , '%S' , SU), '%D', DATE);
    fprintf([' '''  PATH_MAT(34:end-1)  '''\n'])
    end



%%

Sujetos={
 'ASCS_31121976/EEG/08012022_pro'
 'BCUV_02091982/EEG/10122021_pro'
 'C_DG_15111995/EEG/30062021_pro'
 'CACB_27091965/EEG/17122021_pro'
 'CAMO_08111988/EEG/29092021_pro'
 'CARM_14011974/EEG/19032021_pro'
 'CDPR_25031986/EEG/29072021_pro'
 'CEBT_30091989/EEG/03112021_pro'
 'CMCO_17051975/EEG/06102021_pro'
 'EERV_01041955/EEG/25082021_pro'
 'ELFN_09021961/EEG/09092021_pro'
 'FABN_19111984/EEG/19082021_pro'
 'FJAM_10012001/EEG/26112021_pro'
 'GADM_08111983/EEG/26112021_pro'
 'GAQJ_21061987/EEG/16062021_pro'
 'IABF_16091997/EEG/12032021_pro'
 'IABM_03061982/EEG/19082021_pro'
 'JAPV_03111995/EEG/23062021_pro'
 'JAZV_27081991/EEG/05082021_pro'
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
 'VAPH_11101985/EEG/17122021_pro'   
}';




%%
cfg=[];
cfg.subject{1}  =  Sujetos;
cfg.groupname =  {'DATA'};

% %   conditionname = {'task', 'rest', ...}
cfg.    comp 	=[1 ];  	% INDEX OF THE CONDITION TO COMPARED
cfg.    group    =[1 ];  % INDEX OF THE CONDITION TO COMPARED
%     to compare two contrast
%       comp = [ c1 c2 c3 c4]
%       group = [g1 g1 g2 g2]
%       matdif   = [1 -1 2 -2]
%       matdif_transform = 'log', 'none', 'log10'
%       
 cfg.stata='glm';  
  %cfg.glm_matrix = [ones(size(Sujetos))']; ANOSMIA HOSPITAL
  cfg.glm_matrix = [ones(size(Sujetos))'];
  cfg.RegressorOI=1;
 
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
cfg.    range = [5 35];
%cfg.time = [0.5  1.5];
%    permute_name = { 'Real', 'Permut1', 'Permutex', ...}
cfg.matname  = 'LAN';       %                          special carater:
%   .................                                     %S subjectname
cfg.filename  = '/Volumes/DBNC_03/neuroCOVID/DATA/%S/LAN_BL_PUPIL_MODEL.mat';                                     %G groupname 
                                                         %C conditionname  
                                                         %R permute indicator
 cfg.mcp_fast = 0.2;    
% cfg.tempfilename='_model_EYE_rep1' 
% 

GLAN = timefreq_stata([],cfg);



save MODEL_EYE_BL_glm_new_hosp

%

timefreq_plot(GLAN)

