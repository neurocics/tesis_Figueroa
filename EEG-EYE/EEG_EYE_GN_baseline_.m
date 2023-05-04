
clear
clc
% OK 
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


% no data 7 mala calidad
% 'CJGV_24021957''20102021'

SU = SUs{1}  ;%;

DA = '';

LUGAR = '/Volumes/DBNC_03/neuroCOVID/DATA/';


     if isempty(DA)
        [paso paso2] = unix(['ls  -d '  LUGAR SU '/EEG/*_pro' ]);
         DA = paso2(end-12:end-5);
     end





eeg_file = [ LUGAR SU '/EEG/' SU '_' DA  '_GN.eeg'  ];
if ~(exist(eeg_file,'file')==2)
eeg_file = [LUGAR SU '/EEG/OAPC_01041966_23062021.eeg'];
end
log_file = [ LUGAR SU '/LOG/' SU '_' DA  '_GN-gonogo_lum_sin_cali.log'];%
if ~(exist(log_file,'file')==2)
log_file = [ LUGAR SU '/LOG/' SU '_' DA  '_GN-gonogo_lum.log'];%    
end
EYE_POS =  [ LUGAR SU '/LOG/' SU '_' DA  '_GN_EYE_Pos_GNG_file' ];
EYE_PUPIL =[ LUGAR SU '/LOG/' SU '_' DA  '_GN_EYE_Pupil_GNG_file.txt'];

% filtro para pupila
LowPAS = designfilt('lowpassiir', 'PassbandFrequency', 5, 'StopbandFrequency', 8, 'PassbandRipple', 1, 'StopbandAttenuation', 60, 'SampleRate', 1000);

gap=50;
disp(['loading .... >>> ' SU ' <<< '])
load([ LUGAR SU '/EEG/'  DA '_pro/LAN_GN_EEG_interp.mat'])

% check if there are eye tracking data 
if LAN.nbchan<=64
             %   error('NO eye.tracking data')
end

% save ICA componnete   
ica_select = LAN.ica_select;
ica_sphere = LAN.ica_sphere;
ica_weights= LAN.ica_weights;
ica_del= LAN.ica_del;


% read .eeg DATA 
LAN =lan_read_file(eeg_file);%[ SU '_'   DA '_' TAREA '.eeg'  ]
load chanlocs_65geodesic_HCGSNv10(64)
chanlocs(65) = [];
LAN.chanlocs=chanlocs;


% hight past filter in the continuos data 
   if 1
	    d1 = designfilt('highpassiir','FilterOrder',4, ...
	        'HalfPowerFrequency',0.75,'DesignMethod','butter', 'SampleRate',LAN.srate); % 0.25
	    for t=1:LAN.trials    
	    LAN.data{t} = single(filtfilt(d1,double(LAN.data{t}')))';
        end
        disp(['filter OK ....'])
   end

   % inicio baseline 
   S050 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S 50')));
   % inicio tarea
   S040 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S 40')));
   
   
         
        cfg             = [];
        cfg.filename    = log_file;% '../LOG/OAPC_01041962_23062021_WM-Sternberg_EEG__no_eye_tracking_lum.log';% 
        cfg.type        = 'presentation';
        cfg.est         = [50];
        cfg.resp        = [40] ;
        cfg.rw          = [100000000000] ;      %   (ms)    
   RTlog = rt_read(cfg);
   RTlog1 = rt_del( RTlog,2);
   RTlog2 = rt_del( RTlog,1);  
    
   
        cfg             = [];
        cfg.type        = 'RT';
        cfg.RT          = LAN.RT;
        cfg.est         = [S050 ];
        cfg.resp        = [S040];
        cfg.rw          = [1000000000] ;      %   (ms)
        RT_all = rt_read(cfg);
        
   RTlan = rt_read(cfg);  
   
   if numel(RTlan.est)==1
   if abs(RTlan.rt(1)-RTlog.rt(1))  < abs(RTlan.rt(1)-RTlog.rt(2))  % uar si esta perdiod el segundo baseline : est 50 ()
        RTlan = rt_fixlaten(RTlan, RTlog);
   else %1 % si el priemro el segundo 
       RTlan.est = [RTlog.est(1) RTlan.est(1)];
       RTlan.rt = [RTlog.rt(1) RTlan.rt(1)];
       RTlan.laten = [RTlan.laten(1)-(RTlog.laten(2)-RTlog.laten(1))  RTlan.laten(1)];
       RTlan.latency=RTlan.laten ;
       RTlan.resp = [RTlog.resp(1) RTlan.resp(1)];    
       RTlan = rt_check(RTlan);
   end
   RTlan.est(RTlan.est==50) = S050;
   end
   
   RTlan1 = rt_del( RTlan,2);
   RTlan2 = rt_del( RTlan,1);
   

   
   
        % To Segment the data 
        LAN.RT = RTlan1;
        cfg             = [];
        cfg. source     = 'RT'; 
        cfg. ref        = [S050 ] ;   % code of event mark of epoch (references)
        cfg.epoch       = true;
        cfg.times       = [0 (LAN.RT.rt)/1000+1]; % % time for segmentation 
   LAN1 = lan_latency(LAN, cfg );
   LAN1.RT.OTHER.log_latency=RTlog1.laten;
   LAN1.RT.OTHER.log_rt=RTlog1.rt;
        LAN.RT = RTlan2;
        cfg             = [];
        cfg. source     = 'RT'; 
        cfg. ref        = [S050 ] ;   % code of event mark of epoch (references)
        cfg.epoch       = true;
        cfg.times       = [0 (LAN.RT.rt)/1000+1]; % % time for segmentation 
   LAN2 = lan_latency(LAN, cfg );
   LAN2.RT.OTHER.log_latency=RTlog2.laten;
   LAN2.RT.OTHER.log_rt=RTlog2.rt;
   
   
   
   
   
   
   % read EYE TRACKER
   
   % position 
    filename = [EYE_POS  ];
    EYE_POS_LOG =readtable(filename,...
        'Delimiter','\t','ReadVariableNames',true);
    EYE_POS_LOG(EYE_POS_LOG.pos_time==0,:)=[]; % del no data 
    

    % pupil 
    filename = [EYE_PUPIL  ];
    EYE_PUPIL_LOG =readtable(filename,...
        'Delimiter','\t','ReadVariableNames',true);
    EYE_PUPIL_LOG(EYE_PUPIL_LOG.pos_time==0,:)=[]; % del no data 
    
    
    % re sample eye tracker 1200 to 1000 Hz
    f=6;n=0;men='pre(interpolating/resampling  1200Hz to 1000Hz)';
    bar_wait(n,f,men);n=n+1;
    [ POS_X_1000 time_x]= resample(EYE_POS_LOG.pos_x,EYE_POS_LOG.pos_time,1);bar_wait(n,f,men);n=n+1;
    [ POS_Y_1000 ]= resample(EYE_POS_LOG.pos_y,EYE_POS_LOG.pos_time,1);bar_wait(n,f,men);n=n+1;
    [ PUPIL_X_1000 time_pupil]= resample(EYE_PUPIL_LOG.pupil_x ,EYE_PUPIL_LOG.pos_time,1);bar_wait(n,f,men);n=n+1;
    [ PUPIL_Y_1000 ]= resample(EYE_PUPIL_LOG.pupil_y ,EYE_PUPIL_LOG.pos_time,1);bar_wait(n,f,men);n=n+1;

    PUPIL_X_1000 = filtfilt(LowPAS,PUPIL_X_1000);
    PUPIL_Y_1000 = filtfilt(LowPAS,PUPIL_Y_1000);
  
   % convert no measure in NaN in .log data 
    data_inter =  ~ismember(time_x,EYE_POS_LOG.pos_time);
    POS_X_1000_nan = POS_X_1000;
    POS_Y_1000_nan = POS_Y_1000;
    POS_X_1000_nan(data_inter) = NaN;
    POS_Y_1000_nan(data_inter) = NaN;
    
    data_inter_p =  ~ismember(time_pupil,EYE_PUPIL_LOG.pos_time);
    PUPIL_X_1000_nan = PUPIL_X_1000;
    PUPIL_Y_1000_nan = PUPIL_Y_1000;    
    PUPIL_X_1000_nan(data_inter) = NaN;
    PUPIL_Y_1000_nan(data_inter) = NaN;
    
NL = bwlabel(data_inter_p);
for i = 1:max(NL)
    if sum(NL==i) >100
        ini = find(NL==i,1,'first')-gap;%disp(num2str(ini))
        ini = max(ini,1);
        fin = find(NL==i,1,'last')+gap;%disp(num2str(fin))
        fin= min(fin,numel(data_inter_p));
        disp(['BLINK >> '  num2str(sum(NL==i)) ])
        data_inter_p((ini):(fin)) = 1;
    end
end 

  PUPIL_X_1000(data_inter_p)=[];
  PUPIL_Y_1000(data_inter_p)=[];
  time_pupil(data_inter_p)=[];
  
  [ PUPIL_X_1000 ]= resample(PUPIL_X_1000 ,time_pupil,1);bar_wait(n,f,men);n=n+1;
  [ PUPIL_Y_1000 ]= resample(PUPIL_Y_1000,time_pupil,1);bar_wait(n,f,men);n=n+1;
 
   
   
    % Add EYE-TRACKER 

    %syn_log  = time_x(1);

    % fisrt segment  
    LAN=merge_lan(LAN1,LAN2);
    tt = timelan(LAN);
    nchan = 64;%LAN.nbchan;
    for t = 1:length(LAN.data)
        zero_log = find(time_x==fix(LAN.RT.OTHER.log_latency(t)));
        X_log = (zero_log+(LAN.time(t,1)*LAN.srate)):(zero_log+(LAN.time(t,2)*LAN.srate));
        X_log =  X_log(1:end-1);
        red=3;
        amp=30;
        LAN.data{t}(nchan+1,:) = POS_X_1000_nan(X_log)/red;
        LAN.data{t}(nchan+2,:) = POS_X_1000(X_log)/red;
        
        LAN.data{t}(nchan+3,:) = POS_Y_1000_nan(X_log)/red;
        LAN.data{t}(nchan+4,:) = POS_Y_1000(X_log)/red;
        
        LAN.data{t}(nchan+9,:) = PUPIL_Y_1000_nan(X_log)*amp;
        LAN.data{t}(nchan+10,:) = PUPIL_Y_1000(X_log)*amp;
        
        LAN.data{t}(nchan+11,:) = PUPIL_X_1000_nan(X_log)*amp;
        LAN.data{t}(nchan+12,:) = PUPIL_X_1000(X_log)*amp;
    end
    LAN.chanlocs(nchan+1).labels = 'X_nan' ; LAN.chanlocs(nchan+1).type = 'EYE' ;
    LAN.chanlocs(nchan+2).labels = 'X_log' ; LAN.chanlocs(nchan+2).type = 'EYE' ;
    LAN.chanlocs(nchan+3).labels = 'Y_nan' ; LAN.chanlocs(nchan+3).type = 'EYE' ;
    LAN.chanlocs(nchan+4).labels = 'Y_log' ; LAN.chanlocs(nchan+4).type = 'EYE' ;
    LAN.chanlocs(nchan+5).labels = 'XL_tobi' ; LAN.chanlocs(nchan+5).type = 'EYE' ;
    LAN.chanlocs(nchan+6).labels = 'YL_log' ; LAN.chanlocs(nchan+6).type = 'EYE' ;
    LAN.chanlocs(nchan+7).labels = 'XR_tobi' ; LAN.chanlocs(nchan+7).type = 'EYE' ;
    LAN.chanlocs(nchan+8).labels = 'YR_log' ; LAN.chanlocs(nchan+8).type = 'EYE' ;
    LAN.chanlocs(nchan+9).labels = 'P_R_nan' ; LAN.chanlocs(nchan+9).type = 'EYE' ;
    LAN.chanlocs(nchan+10).labels = 'P_R_log' ; LAN.chanlocs(nchan+10).type = 'EYE' ; % Y is Rigth 
    LAN.chanlocs(nchan+11).labels = 'P_L_nan' ; LAN.chanlocs(nchan+11).type = 'EYE' ;
    LAN.chanlocs(nchan+12).labels = 'P_L_log' ; LAN.chanlocs(nchan+12).type = 'EYE' ; % X is left 
     LAN = lan_check(LAN);

%
close all
% mark bad segment 
prepro_plot(LAN)


%% segemntar y preprosesar %   
LAN.ica_select = ica_select;
LAN.ica_sphere = ica_sphere;
LAN.ica_weights= ica_weights;
%relet  ICA
LAN = lan_rm_chan(LAN, ica_del,'ica');

LAN = lan_segment_selected(LAN,2,1);

elec = 1:58;
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
prepro_plot(LAN) % select trials / e



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

%


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


%

save([ LUGAR SU '/EEG/'  DA '_pro/LAN_GN_EEG_EYE_Baseline.mat'],'LAN')
disp(['DONE ' SU ' OK'])
clear