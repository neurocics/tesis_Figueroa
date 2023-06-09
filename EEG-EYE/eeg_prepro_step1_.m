
% This script DO:   
%               1. Read .egg data fron BRAINAMP 
%               2. Read .log dat fron PRESENTION 
%               3. Synchrionize 1 and 2 -  Recover lost marker 
%               4. Segmente LAN (.eeg data) 
%               5. Synchronized EYE TRACKER data 
%               6. Added to LAN structure and save it 
%                  
% 
% dependency:
% LAN toolbox   <*LAN)<] v.1.9.9
% 
%  15 Sept 2021
%  Pablo Billeke 

clc
tic %  
clear

% SET variables per subjet 

% names of subejto to build  file and directories
GR = 'controles';
SU = 'CFPL_08021984';
DATA = '02052023';
TAREA = 'RL';

% cases:  'WM' for working memory , 180 TRIALS
%         'RL' for reversal learning 
%         'GN' for Go no go 





% eyetracking file 
EYE_TOBI =[SU '_' DATA  '_' TAREA '.txt'];% 'AGBB_26121972.txt'%
RI=1;

% eeg file 
eeg_file =[ SU '_'   DATA '_' TAREA '.eeg'  ];% 'GAQJ_20061987_16062021_RL.eeg'%'LAAB_29101982_02062021_RL.eeg';%

% recording index from the EYE tracker file




% file per task 
% CHECK the specific speeling for each subject 
switch TAREA
     case 'GN' 
        % 
        %LOGP = '_GN-gonogo_no_EYE_lum.log'; EYE=0;LOG=0;TOBI=0;% for no EYE-TRACKING data
        LOGP = '_GN-gonogo_lum_sin_cali.log';EYE=1;LOG=1;TOBI=0; %
        EYE_POS ='_GN_EYE_Pos_GNG_file.txt';
        EYE_PUPIL ='_GN_EYE_Pupil_GNG_file.txt';
        SEG_TIME = [-1.5 1.5];  %  
    case 'WM'
        % % for no EYE-TRACKING data
        %LOGP = '-Sternberg_EEG__no_eye_tracking_lum.log'; EYE=0;LOG=1;TOBI=0;
        LOGP = '-Sternberg_EEG__eye_tracking_lum_sin_cali.log';EYE=1;LOG=1;TOBI=0; %
        EYE_POS ='_WM_EYE_Pos_WM_file.txt';
        EYE_PUPIL ='_WM_EYE_Pupil_WM_file.txt';
        SEG_TIME = [-1 5.5];  % time for segemetation 0: array, 1.8: end encoiding, 4: end mantention/probe  
    case 'RL'
        % 
        %LOGP = '_RL-RvL_no_eye.log';EYE=0;LOG=1;TOBI=0; % for no EYE-TRACKING data
        LOGP = '_RL-RvL_eye_sin_cali.log';EYE=1;LOG=1;TOBI=0;%
        EYE_POS ='_RL_EYE_Pos_RvL_file.txt';
        EYE_PUPIL ='_RL_EYE_Pupil_RvL_file.txt';
        OUT = '_RL_output_RvL_file.txt';
        SEG_TIME = [-1.5 4];  % the last in the minum, ...
        % the scrip adat to segment all the decision time fot the 95% of
        % the trial, thus each subject
        % have different lenght... in subsequent analisis it must be fixed  
        % in order to get a constant length of the trials for all subject 
end






% directories %S = subject; %G=groups (time of evaluation)
%PATH_LOG =  '/Volumes/DB_neuroCICS_02/datos_investigadores/neuroCOVID/%G/%S/LOG/';
%PATH_LOG = '/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/%S/LOG/';
%PATH_TOBI =  '/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/%S/EYE/';
%PATH_D =  '/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/%S/EEG/';
%PATH_LOG = '/Volumes/LaCie/TESIS_FIGUEROA/DATA/%G/%S/LOG/';
%PATH_TOBI = '/Volumes/LaCie/TESIS_FIGUEROA/DATA/%G/%S/EYE/';
%PATH_D =  '/Volumes/LaCie/TESIS_FIGUEROA/DATA/%G/%S/EEG/';
PATH_LOG = '/Volumes/Alehermosa/TESIS/%G/%S/LOG/';
PATH_TOBI = '/Volumes/Alehermosa/TESIS/%G/%S/EYE/';
PATH_D =  '/Volumes/Alehermosa/TESIS/%G/%S/EEG/';





% create real PATH for the files 
PATH_LOG =strrep( strrep(PATH_LOG,'%G',GR) , '%S' , SU);
PATH_TOBI =strrep( strrep(PATH_TOBI,'%G',GR) , '%S' , SU);
PATH_D = strrep( strrep(PATH_D,'%G',GR) , '%S' , SU);


%% STEPS: 1 to 3  


% create directoy to save proceced DATA
cd(PATH_D);
%mkdir MAT

% read .egg DATA 
LAN =lan_read_file(eeg_file(1:end-4),'BA');%[ SU '_'   DATA '_' TAREA '.eeg'  ]

%LAN =lan_read_file(eeg_file);%[ SU '_'   DATA '_' TAREA '.eeg'  ]
% % fin no label % posible error using FILEIO toolbox
% no_label = ifcellis(LAN.RT.OTHER.names,'_no_label');
% paso = nan(size(LAN.RT.est));
% paso(~no_label) = LAN.RT.laten;
% paso(no_label) = 1;
% LAN.RT.laten = paso;
% LAN.RT.latency = paso;

% hight past filter in the continuos data 
   if 1    % cambia mejor a FIR 
	    d1 = designfilt('highpassiir','FilterOrder',4, ...
	        'HalfPowerFrequency',0.75,'DesignMethod','butter', ...
            'SampleRate',LAN.srate); % 0.25

       % d1 = designfilt('highpassfir','FilterOrder',4,'StopbandFrequency',0.2, ...
       %                'PassbandFrequency',0.75,'SampleRate',LAN.srate);


	    for t=1:LAN.trials    
	    LAN.data{t} = single(filtfilt(d1,double(LAN.data{t}')))';
        end
        disp(['filter OK ....'])
	end


switch TAREA
    case 'GN'
        %-----------------------
    %---%  GO NO GO TASK 
        %-----------------------

        % Stimuli in .egg data 
 
        S010 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S 10')));% GO
        S021 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S 21')));% NO GO
        % S050 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S 50')));
        % % BASE LINE 
        %[S010 S021 S050]
        % responses
        S001 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S  1')));
        S002 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S  2')));
        %[S001 S002]


        % extract to each stim the corresponding probe and response form .eeg
        % and .log of presentation 
        %
        % .eeg
        cfg             = [];
        cfg.type        = 'RT';
        cfg.RT          = LAN.RT;
        cfg.est         = [S001 S002 S010 S021 ];
        cfg.resp        = [200];
        cfg.rw          = [1] ;      %   (ms)
        RT_all = rt_read(cfg);
        
        if LOG
        % .log
        cfg             = [];
        cfg.filename    = ['../LOG/' SU '_' DATA   LOGP ];% '../LOG/OAPC_01041962_23062021_WM-Sternberg_EEG__no_eye_tracking_lum.log';% 
        cfg.type        = 'presentation';
        cfg.est         = [1 2 10 21];
        cfg.resp        = [200] ;
        cfg.rw          = [1] ;      %   (ms)    
        RT_all_log = rt_read(cfg);

        log_code = RT_all_log.est;

        % Change for an equivalente CODE between RT form .log adn RT from .eeg 
        RT_all_log.est(RT_all_log.est==1) = S001;
        RT_all_log.est(RT_all_log.est==2) = S002;

        RT_all_log.est(RT_all_log.est==10) = S010;
        RT_all_log.est(RT_all_log.est==21) = S021;
        end
    case 'WM'
        %-----------------------
    %---%  WORKING MEMEORY TASK 
        %-----------------------

        % Stimuli in .egg data 
        % array stim 
        S012 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S 12')));
        S014 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S 14')));
        S016 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S 16')));
        %[S012 S014 S016]
        % kind of probe
        S020 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S 20')));
        S021 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S 21')));
        S040 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S 40')));
        S041 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S 41')));
        S060 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S 60')));
        S061 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S 61')));
        %[S020 S021 S040 S041 S060 S061]
        % responses
        S001 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S  1')));
        S002 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S  2')));
        %[S001 S002]


        % extract to each stim the corresponding probe and response form .eeg
        % and .log of presentation 
        %
        % .eeg
        cfg             = [];
        cfg.type        = 'RT';
        cfg.RT          = LAN.RT;
        cfg.est         = [S001 S002 S012 S014 S016 S020 S021 S040 S041 S060 S061];
        cfg.resp        = [200];
        cfg.rw          = [1] ;      %   (ms)
        RT_all = rt_read(cfg);
        if LOG
        % .log
        cfg             = [];
        cfg.filename    = ['../LOG/' SU '_' DATA '_' TAREA  LOGP ];% '../LOG/OAPC_01041962_23062021_WM-Sternberg_EEG__no_eye_tracking_lum.log';% 
        cfg.type        = 'presentation';
        cfg.est         = [1 2 12 14 16 20 21 40 41 60 61];
        cfg.resp        = [200] ;
        cfg.rw          = [1] ;      %   (ms)    
        RT_all_log = rt_read(cfg);
        
        log_code = RT_all_log.est;

        % Change for an equivalente CODE between RT form .log and RT from
        % .eeg (equivalencia entre estimulos de presentation y eeg)
        
        RT_all_log.est(RT_all_log.est==1) = S001;
        RT_all_log.est(RT_all_log.est==2) = S002;

        RT_all_log.est(RT_all_log.est==12) = S012;
        RT_all_log.est(RT_all_log.est==14) = S014;
        RT_all_log.est(RT_all_log.est==16) = S016;

        RT_all_log.est(RT_all_log.est==20) = S020;
        RT_all_log.est(RT_all_log.est==21) = S021;
        RT_all_log.est(RT_all_log.est==40) = S040;
        RT_all_log.est(RT_all_log.est==41) = S041;
        RT_all_log.est(RT_all_log.est==60) = S060;
        RT_all_log.est(RT_all_log.est==61) = S061;
        end
    case 'RL'
        %-----------------------
    %---%  REVERSAL LEARNING TASK 
        %-----------------------     
        % Stimuli in .egg data 
        % opstions  
        S010 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S 10')));
        S015 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S 15')));
        %[S010 S015]
        % feedback
        S030 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S 30')));
        S031 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S 31')));
        %[S030 S031]
        % responses
        S001 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S  1')));
        S002 = unique(LAN.RT.est(ifcellis(LAN.RT.OTHER.names , 'S  2')));
        %[S001 S002]

        % a case with no stimuli event     
        if isempty(S010)&&isempty(S015)&&isempty(S030)&&isempty(S031)
           S010=10;S015=15;S030=30;S031=31;    
        end

        % extract to each stim the corresponding probe and response form .eeg
        % and .log of presentation 
        %
        % .eeg
        cfg             = [];
        cfg.type        = 'RT';
        cfg.RT          = LAN.RT;
        cfg.est         = [S001 S002 S010 S015 S030 S031];
        cfg.resp        = [200];
        cfg.rw          = [1] ;      %   (ms)
        RT_all= rt_read(cfg);

        % .log
        cfg             = [];
        cfg.filename    = ['../LOG/' SU '_' DATA   LOGP ];% '../LOG/OAPC_01041962_23062021_WM-Sternberg_EEG__no_eye_tracking_lum.log';% 
        cfg.type        = 'presentation';
        cfg.est         = [1 2 10 15 30 31];
        cfg.resp        = [200] ;
        cfg.rw          = [1] ;      %   (ms)    
        RT_all_log = rt_read(cfg);
%         if LOG
%            RT_all_log =rt_del(RT_all_log,find(RT_all_log.est<3)+1);
%         end
        log_code = RT_all_log.est;
        
                % Change for an equivalente CODE between RT form .log adn RT from .eeg 
        RT_all_log.est(RT_all_log.est==1) = S001;
        RT_all_log.est(RT_all_log.est==2) = S002;
        %----- revisar esto bien 
        try 
        try
        RT_all_log.est(RT_all_log.est==10) = S010;
        catch 
        RT_all_log.est(RT_all_log.est==10) = S010;    
        end
        RT_all_log.est(RT_all_log.est==15) = S015;
        %-------------------------
       
        RT_all_log.est(RT_all_log.est==30) = S030;
        RT_all_log.est(RT_all_log.est==31) = S031;
        end
        
        
end 

if 0
figure
% util para ver que esta mal 
plot(RT_all.est,'r'), hold on, plot(RT_all_log.est,'k')
[RT_all.est(1:16)' RT_all_log.est(1:16)']
end






% % fix missing latency
if LOG
RT_all_log.dw_delta=30;
RT_all_log.up_delta=30;
%

% % %RT_all_fix = rt_fixlaten(RT_all,rt_del(RT_all_log,1:11));
% RT_all=rt_del(RT_all,1:2);
%  RT_all_log = rt_del(RT_all_log,1:11);%

%[ RT_all.est(1:30)'  RT_all_log.est(1:30)']


RT_all_fix = rt_fixlaten(RT_all,RT_all_log);
% trata de corregir entre informacion enviada de EEG y presentation
% RT_all_fix = rt_del(RT_all_fix,1:2);

% lost trials in the eeg recording for error in the recording 
n_lost_eeg_trial=0;


% CHECK if it is all OK !!!
% in this point can be an error if the first stimuli is lost 
% try this 

%P = rt_del(RT_all,1)

%[ RT_all_fix.laten(RT_all_fix.est<3)' P.laten(P.est<3)'  ]
%[ RT_all.est(1:20)' P.laten(P.est<3)'  ]
% RT_all_fix = rt_fixlaten(rt_del(RT_all,1),rt_del(RT_all_log,1));


%  lost_stim = 1:68;
%  % rt_fixlaten(RT_all,rt_del(RT_all_log,lost_stim));
% % 
%  RT_all_log = rt_del(RT_all_log,lost_stim);
%  RT_all_fix = rt_fixlaten(RT_all,RT_all_log); %
%  log_code = RT_all_log.est;
end


% Save the time of the stimuli from .log in order to sychronize EYE-TRACKER
% data 
if LOG
RT_all_fix.OTHER.log_latency = RT_all_log.latency;
RT_all_fix.OTHER.log_code = log_code;
% rempace the RT from the fixed one
LAN.RT = RT_all_fix;
end



%% STEP: 4  use the fixed RT to segemnet .eeg data 
switch TAREA
    case 'GN'
        %-----------------------
    %---%  GO NO GO  TASK 
        %-----------------------
        
        % go no-go stimuli
        % .eeg
        cfg             = [];
        cfg.type        = 'RT';
        cfg.RT          = LAN.RT;
        cfg.est         = [S010 S021];
        cfg.resp        = [S001 S002] ;
        cfg.rw          = [10000] ;      %   (ms)
        RT_stim = rt_read(cfg);
        % .log
        if LOG
        cfg             = [];
        cfg.filename    = ['../LOG/' SU '_' DATA   LOGP ];
        cfg.type        = 'presentation';
        cfg.est         = [10 21];
        cfg.resp        = [1 2] ;
        cfg.rw          = [10000] ;      %   (ms)    %   (ms)
        RT_stim_log = rt_read(cfg);
        end
        % To Segment the data 
        cfg             = [];
        cfg. source     = 'RT'; 
        cfg. ref        = [S010 S021] ;   % code of event mark of epoch (references)
        cfg.resp        = [S001 S002] ;
        cfg.epoch       = true;
        cfg.times       = SEG_TIME; % % time for segmentation 
        LAN = lan_latency(LAN, cfg );
            if LOG
            % check vector lenght 
            if  ~(  length(RT_stim.est)==length(RT_stim_log.est)...
                    )
               error([' Error in length of  RT vectors !!!']);
            else
                fprintf(['\n\t\tEverything is OK .... <*LAN)< \n\n'])
            end
            end

        % add to each stim the corresponding probe and response
        %LAN.RT.rt = RT_stim.rt;
        %LAN.RT.resp = RT_stim.resp;
        if LOG
        LAN.RT.OTHER.log_RT = RT_stim_log.rt;
        LAN.RT.OTHER.log_resp = RT_stim_log.resp;
        end
    case 'WM'
        %-----------------------
    %---%  WORKING MEMEORY TASK 
        %-----------------------
        
        % the kind of probe for each stimulus
        % .eeg
        cfg             = [];
        cfg.type        = 'RT';
        cfg.RT          = LAN.RT;
        cfg.est         = [S012 S014 S016];
        cfg.resp        = [S020 S021 S040 S041 S060 S061] ;
        cfg.rw          = [10000] ;      %   (ms) hasta cuando busca una respuesta del sujeto
        RT_probe = rt_read(cfg);
        % .log
        if LOG
        cfg             = [];
        cfg.filename    = ['../LOG/' SU '_' DATA '_' TAREA  LOGP ];
        cfg.type        = 'presentation';
        cfg.est         = [12 14 16 ];
        cfg.resp        = [20 21 40 41 60 61] ;
        cfg.rw          = [10000] ;      %   (ms)    %   (ms)
        RT_probe_log = rt_read(cfg);
        if n_lost_eeg_trial
            RT_probe_log = rt_del(RT_probe_log,1:n_lost_eeg_trial) 
        end
        end

        % the response for each stimulus
        cfg             = [];
        cfg.type        = 'RT';
        cfg.RT          = LAN.RT;
        cfg.est         = [S012 S014 S016];
        cfg.resp        = [S001 S002] ;
        cfg.rw          = [10000] ;      %   (ms)
        RT_responses = rt_read(cfg);

        % To Segment the data 
        cfg             = [];
        cfg. source     = 'RT'; 
        cfg. ref        = [S012 S014 S016] ;   % code of event mark of epoch (references)
        cfg.epoch       = true;
        cfg.times       = SEG_TIME; % % time for segmentation 
        LAN = lan_latency(LAN, cfg );
            if LOG
            % check vector lenght 
            if  ~(  length(RT_probe.est)==length(RT_responses.est) &&...
                    length(RT_probe.est)==length(RT_probe_log.resp) &&...
                    length(RT_probe.est)==length(LAN.RT.est)...
                    )
               error([' Error in length of  RT vectors !!!']);
            else
                fprintf(['\n\t\tEverything is OK .... <*LAN)< \n\n'])
            end
            end
        % add to each stim the corresponding probe and response
        LAN.RT.rt = RT_responses.rt;
        LAN.RT.resp = RT_responses.resp;
        LAN.RT.OTHER.probe = RT_probe.resp;
        if LOG
        LAN.RT.OTHER.log_probe = RT_probe_log.resp;
        end
    case 'RL'
        %-----------------------
    %---%  REVERSAL LEARNING TASK 
        %-----------------------
        
        % the response  for each stimulus
        % .eeg
        cfg             = [];
        cfg.type        = 'RT';
        cfg.RT          = rt_del( LAN.RT, find(LAN.RT.good==0));
        cfg.est         = [S010 S015 ];
        cfg.resp        = [S001 S002] ;
        cfg.rw          = [200000] ;      %   (ms)
        RT_responses = rt_read(cfg);
        %si = [find(RT_responses.resp>0)-1];
        
        
        % delete stimuli thata there are no option presentation (e.i., the begining of the waiting time )
         RT_responses = rt_del(RT_responses , 2:2:(numel(RT_responses.est)));
        % .log
        cfg             = [];
        cfg.filename    = ['../LOG/' SU '_' DATA   LOGP ];
        cfg.type        = 'presentation';
        cfg.est         = [10 15 ];
        cfg.resp        = [1 2] ;
        cfg.rw          = [20000] ;      %   (ms)    %   (ms)
        RT_responses_log = rt_read(cfg);
        % delete stimuli thata there are no option presentation (e.i., the begining of the waiting time )        
        RT_responses_log = rt_del(RT_responses_log , 2:2:(numel(RT_responses_log.est)));
        
        % out file with trial information 

        OUT_FILE =readtable(['../LOG/' SU '_' DATA   OUT],...
         'Delimiter','\t','ReadVariableNames',true);   

     
        %------------------------------
        % lost the fisrt n trials in the .eeg

        if n_lost_eeg_trial
        OUT_FILE(1:n_lost_eeg_trial,:) =[]; 
        RT_responses_log = rt_del(RT_responses_log,1:n_lost_eeg_trial);
        end

        %-------------------------------
        
        % check vector lenght 
            if  ~(...
                    length(RT_responses.est)==length(RT_responses_log.est) && ...
                    length(RT_responses_log.est)==length(OUT_FILE.nt)...
                )
               error([' Error in length of  RT vectors for decision  !!!']);
            else
                fprintf(['\n\t\t decision stimuli are  OK ....  <*LAN)< \n\n'])
            end
        % add to each stim the corresponding log and out-file information 
        RT_responses.OTHER.log_latency = RT_responses_log.latency;
        RT_responses.OTHER.out_latency = OUT_FILE.ttp;    
        RT_responses.OTHER.out_rt = OUT_FILE.ttr - OUT_FILE.ttp;  
        RT_responses.OTHER.out_Bd1 = OUT_FILE.Bd1;	
        RT_responses.OTHER.out_Bd2 = OUT_FILE.Bd2;	        
        RT_responses.OTHER.out_buena = OUT_FILE.buena;	       
        RT_responses.OTHER.out_resp = OUT_FILE.Resp;
        RT_responses.OTHER.out_F = OUT_FILE.F;
        RT_responses.OTHER.out_prob = OUT_FILE.prob;
        RT_responses.OTHER.out_nt = OUT_FILE.nt;
        
        % the feedback  for each stimulus
        % .eeg
        cfg             = [];
        cfg.type        = 'RT';
        cfg.RT          = LAN.RT;
        cfg.est         = [S030 S031];
        cfg.resp        = [S010 S015] ; % non resposnes , only the posible following stimuli 
        cfg.rw          = [20000] ;      %   (ms)
        RT_feedback = rt_read(cfg);
        % .log
        cfg             = [];
        cfg.filename    = ['../LOG/' SU '_' DATA   LOGP ];
        cfg.type        = 'presentation';
        cfg.est         = [30 31 ];
        cfg.resp        = [10 15] ; % non resposnes , only the posible following stimuli 
        cfg.rw          = [20000] ;      %   (ms)    %   (ms)
        RT_feedback_log = rt_read(cfg);
        
        %------------------------------
        % lost the fisrt n trials in the .eeg
        % % n_lost_eeg_trial=1;
        if n_lost_eeg_trial
        RT_feedback_log = rt_del(RT_feedback_log,1:n_lost_eeg_trial);
        % for inclomple first trial!
        if length(RT_feedback.est) == length(RT_feedback_log.est) +1;
            RT_feedback = rt_del(RT_feedback,1:n_lost_eeg_trial);
        end
        end
        %-------------------------------
        
        
        
        % check vector lenght 
            if  ~(...
                    length(RT_feedback.est)==length(RT_feedback_log.est) && ...
                    length(RT_feedback_log.est)==length(OUT_FILE.nt)...
                )
               error([' Error in length of  RT vectors for feedback  !!!']);
            else
                fprintf(['\n\t\t feedback stimuli are  OK ....  <*LAN)< \n\n'])
            end
                % add to each stim the corresponding log and out-file information 
        RT_feedback.OTHER.log_latency = RT_feedback_log.latency;
        RT_feedback.OTHER.out_latency = OUT_FILE.ttF;    
        RT_feedback.OTHER.out_rt = OUT_FILE.ttr - OUT_FILE.ttp;  
        RT_feedback.OTHER.out_Bd1 = OUT_FILE.Bd1;	
        RT_feedback.OTHER.out_Bd2 = OUT_FILE.Bd2;	        
        RT_feedback.OTHER.out_buena = OUT_FILE.buena;	       
        RT_feedback.OTHER.out_resp = OUT_FILE.Resp;
        RT_feedback.OTHER.out_F = OUT_FILE.F;
        RT_feedback.OTHER.out_prob = OUT_FILE.prob;
        RT_feedback.OTHER.out_nt = OUT_FILE.nt;
        % Merge RT and sorting  stimli 
        LAN.RT = rt_merge(RT_responses, RT_feedback,1);
       
        % % adat to segmente all the decision time plus 0.5s (limit 6 secong )
         max_time = max(RT_responses.rt); 
         per = sort(RT_responses.rt);
         disp(['max time      : '  num2str(max(RT_responses.rt))]);
         disp(['99 time      : '  num2str(per(fix(numel(per)*.99)))]);
         disp(['97 time      : '  num2str(per(fix(numel(per)*.97)))]);
         disp(['95 time      : '  num2str(per(fix(numel(per)*.95)))]);
         disp(['90 time      : '  num2str(per(fix(numel(per)*.90)))]);
         disp(['80 time      : '  num2str(per(fix(numel(per)*.80)))]);
         disp(['70 time      : '  num2str(per(fix(numel(per)*.70)))]);
        if SEG_TIME(2)<(0.5 + per(fix(numel(per)*.95))/1000)
           SEG_TIME(2)=(0.5 + per(fix(numel(per)*.95))/1000);
        end
         
        % To Segment the data 
        cfg             = [];
        cfg. source     = 'RT'; 
        cfg. ref        = [S010 S015 S030 S031 ] ;   % code of event mark of epoch (references)
        cfg.epoch       = true;
        cfg.times       = SEG_TIME; % % time for segmentation 
        LAN = lan_latency(LAN, cfg );


 
end


%% STEP: 5
% Extracted and synchronized de EYE tracker recordings 

if EYE

    % read EYE tracker data readout form presentation 
    % --> necesary for the synchronization 

    % position 
    filename = [PATH_LOG SU '_' DATA  EYE_POS  ];
    EYE_POS_LOG =readtable(filename,...
        'Delimiter','\t','ReadVariableNames',true);
    EYE_POS_LOG(EYE_POS_LOG.pos_time==0,:)=[]; % del no data 
    

    % pupil 
    filename = [PATH_LOG SU '_' DATA  EYE_PUPIL  ];
    EYE_PUPIL_LOG =readtable(filename,...
        'Delimiter','\t','ReadVariableNames',true);
    EYE_PUPIL_LOG(EYE_PUPIL_LOG.pos_time==0,:)=[]; % del no data 
    
    
    % fix laten for each pauses in .log 
    % MEJORAR ESTO, POR QUE LA INTERPOLACION FALLA 
    if isfield(RT_all_log,'PAUSES') && ~isempty(RT_all_log.PAUSES) 
      
    dif_t =[0 ; diff(EYE_POS_LOG.pos_time)];
    lag_tim = find( dif_t <0)' ;
    for il =  1:numel(lag_tim)
       lag_tim_i(il) = find_approx(dif_t,abs(dif_t(lag_tim(il)))   ); 
    end
        for il =  numel(lag_tim):-1:1
        EYE_POS_LOG(lag_tim_i(il):lag_tim(il),:) =[]  ; 
        end
    
    
    dif_t =[0 ; diff(EYE_PUPIL_LOG.pos_time)];
    lag_tim = find( dif_t <0)' ;
    for il =  1:numel(lag_tim)
       lag_tim_i(il) = find_approx(dif_t,abs(dif_t(lag_tim(il)))   ); 
    end
        for il =  numel(lag_tim):-1:1
        EYE_PUPIL_LOG(lag_tim_i(il):lag_tim(il),:) =[]  ; 
        end
   
    
    
    
%     plot(EYE_POS_LOG.pos_time)   , hold on 
%     for il =  1:numel(lag_tim)
%     line([lag_tim_i(il) lag_tim_i(il)], [min(EYE_POS_LOG.pos_time)  max(EYE_POS_LOG.pos_time)] );
%     line([lag_tim(il) lag_tim(il)], [min(EYE_POS_LOG.pos_time)  max(EYE_POS_LOG.pos_time)] );
%     end


        
      for ip = size(RT_all_log.PAUSES,1):-1:1
          EYE_POS_LOG.pos_time(EYE_POS_LOG.pos_time>(RT_all_log.PAUSES(ip,2))) = EYE_POS_LOG.pos_time(EYE_POS_LOG.pos_time>(RT_all_log.PAUSES(ip,2))) + RT_all_log.PAUSES(ip,3);
          EYE_PUPIL_LOG.pos_time(EYE_PUPIL_LOG.pos_time>(RT_all_log.PAUSES(ip,2))) = EYE_PUPIL_LOG.pos_time(EYE_PUPIL_LOG.pos_time>(RT_all_log.PAUSES(ip,2))) + RT_all_log.PAUSES(ip,3);
      end
    end    
    
  
    
    
    if TOBI
        
               
    filename = [PATH_TOBI  EYE_TOBI  ];
    
    
      if exist(filename,'file')~=2
         filename2 = [ filename(1:end-3)  'tsv'];
         if exist(filename2,'file')==2
             filename_f = strrep(filename,' ','\ ');filename2_f = strrep(filename2,' ','\ ');
             r = unix(['cp ' filename2_f ' ' filename_f ' ']);
             if r == 0
                 unix(['rm ' filename2_f ' '])
             end
         end
      end
    
    %Data=importdata('DM.txt');
    EYE_ALL_TOBI =readtable(filename,...
        'Delimiter','\t','ReadVariableNames',true);
    
    %EYE_ALL_TOBI(find(diff(EYE_ALL_TOBI.RecordingTimestamp)==0),:) = [];

    % read EYE tracker data readout form tobii prolab 
    % --> necesary for aditional information 
      
    %     % % evalaute the oreder osf the task // RECOMENDE: one file per task !!!
    %       figure, plot(EYE_ALL_TOBI.RecordingTimestamp), hold on
    %       plot(EYE_ALL_TOBI.ComputerTimestamp)
    %       figure,plot(EYE_ALL_TOBI.GazePointX)
    %      
    %        figure,plot(EYE_POS_LOG.pos_x)
    %     % evalaute for each task 
    % recording index
    % RI=2;  % setting in the HEADER
    
    [index] = find(diff(EYE_ALL_TOBI.RecordingTimestamp)<0);
    index= [1 index' numel(EYE_ALL_TOBI.RecordingTimestamp)];
    index = (index(RI)+1):(index(RI+1)-1);
    
    %------------------------------------------------------------------
    %cut extrem 
    %    when the extrem are too noises, this cut improbe ssynchronization
    %cut_end = 1.1e6;% 
    cut_end = numel(index);%
    %cut_begin = 0.65e6;%
    cut_begin = 1;%
    index = index(cut_begin:cut_end);
    %------------------------------------------------------------------
    end   
    
    % use the following if there is only one task in the file: it is
    % RECOMENDED 
    % index =1:numel(EYE_ALL_TOBI.RecordingTimestamp)

    % expant time in GAPS with missing data 
    %  linear interpolation/resample  to 1000HZ (to be compatible with .eeg data)
    %  this part it is low, and may be no necesary in this way.
   
    
    % for waiting bar 
    f=20;n=0;men='pre(interpolating/resampling  1200Hz to 1000Hz)';
    bar_wait(n,f,men);n=n+1;
   
  
    
    % .log
   % if TOBI
    [ POS_X_1000 time_x]= resample(EYE_POS_LOG.pos_x,EYE_POS_LOG.pos_time,1);bar_wait(n,f,men);n=n+1;
    [ POS_Y_1000 ]= resample(EYE_POS_LOG.pos_y,EYE_POS_LOG.pos_time,1);bar_wait(n,f,men);n=n+1;
    [ PUPIL_X_1000 time_pupil]= resample(EYE_PUPIL_LOG.pupil_x ,EYE_PUPIL_LOG.pos_time,1);bar_wait(n,f,men);n=n+1;
    [ PUPIL_Y_1000 ]= resample(EYE_PUPIL_LOG.pupil_y ,EYE_PUPIL_LOG.pos_time,1);bar_wait(n,f,men);n=n+1;
%     else
%      i_pos =[  diff(EYE_POS_LOG.pos_time)==0 ; 0] | [  abs(diff(EYE_POS_LOG.pos_time))>5000 ; 0]; 
%      POS_X_1000 = EYE_POS_LOG.pos_x(~i_pos) ;
%      POS_Y_1000 = EYE_POS_LOG.pos_x(~i_pos) ;
%      time_x = EYE_POS_LOG.pos_time(~i_pos) ;
%      
%      i_pos =[  diff(EYE_PUPIL_LOG.pos_time)==0 ; 0] | [  abs(diff(EYE_PUPIL_LOG.pos_time))>5000 ; 0]; 
%      PUPIL_X_1000 = EYE_PUPIL_LOG.pupil_x(~i_pos) ;
%      PUPIL_Y_1000 = EYE_PUPIL_LOG.pupil_y(~i_pos) ;
%      time_pupil = EYE_PUPIL_LOG.pos_time(~i_pos) ; 
%      
%     end

    % .tsv Tobii
    % fix precision  in the time recording in tobii file 
    % this is the logic (and faster) way but do not work fine in the sinchronization 
    % fixed_tobii_time = linspace(EYE_ALL_TOBI.RecordingTimestamp(index(1)),EYE_ALL_TOBI.RecordingTimestamp(index(end)),numel(index));
   if TOBI
    fixed_tobii_time =EYE_ALL_TOBI.RecordingTimestamp(index);
    %[ POS_X_1000_tobi time_x_tobi]= resample(EYE_ALL_TOBI.GazePointX(index),fixed_tobii_time,1);bar_wait(n,f,men);n=n+1;
    %[ POS_Y_1000_tobi ]= resample(EYE_ALL_TOBI.GazePointY(index),fixed_tobii_time,1);bar_wait(n,f,men);n=n+1;
    [ POS_XL_1000_tobi ]= resample(EYE_ALL_TOBI.GazePointLeftX(index),fixed_tobii_time,1);bar_wait(n,f,men);n=n+1;
    [ POS_YL_1000_tobi ]= resample(EYE_ALL_TOBI.GazePointLeftY(index),fixed_tobii_time,1);bar_wait(n,f,men);n=n+1;
    [ POS_XR_1000_tobi ]= resample(EYE_ALL_TOBI.GazePointRightX(index),fixed_tobii_time,1);bar_wait(n,f,men);n=n+1;
    [ POS_YR_1000_tobi ]= resample(EYE_ALL_TOBI.GazePointRightY(index),fixed_tobii_time,1);bar_wait(n,f,men);n=n+1;

    [ Gdir_XL_1000_tobi ]= resample(EYE_ALL_TOBI.GazeDirectionLeftX(index),fixed_tobii_time,1);bar_wait(n,f,men);n=n+1;
    [ Gdir_YL_1000_tobi ]= resample(EYE_ALL_TOBI.GazeDirectionLeftY(index),fixed_tobii_time,1);bar_wait(n,f,men);n=n+1;
    [ Gdir_ZL_1000_tobi ]= resample(EYE_ALL_TOBI.GazeDirectionLeftZ(index),fixed_tobii_time,1);bar_wait(n,f,men);n=n+1;

    [ Gdir_XR_1000_tobi ]= resample(EYE_ALL_TOBI.GazeDirectionRightX(index),fixed_tobii_time,1);bar_wait(n,f,men);n=n+1;
    [ Gdir_YR_1000_tobi ]= resample(EYE_ALL_TOBI.GazeDirectionRightY(index),fixed_tobii_time,1);bar_wait(n,f,men);n=n+1;
    [ Gdir_ZR_1000_tobi ]= resample(EYE_ALL_TOBI.GazeDirectionRightZ(index),fixed_tobii_time,1);bar_wait(n,f,men);n=n+1;

    [ PUPIL_L_1000_tobi time_pupil_tobi]= resample(EYE_ALL_TOBI.PupilDiameterLeft(index),fixed_tobii_time,1);bar_wait(n,f,men);n=n+1;
    [ PUPIL_R_1000_tobi]= resample(EYE_ALL_TOBI.PupilDiameterRight(index),fixed_tobii_time,1);bar_wait(n,f,men);n=n+1;

    VALID_L = ifcellis(EYE_ALL_TOBI.ValidityLeft,'Valid');
    VALID_R = ifcellis(EYE_ALL_TOBI.ValidityRight,'Valid');
    [ VALID_L_1000_tobi ]= resample(double(VALID_L(index)),fixed_tobii_time,1);bar_wait(n,f,men);n=n+1;
    [ VALID_R_1000_tobi ]= resample(double(VALID_R(index)),fixed_tobii_time,1);bar_wait(n,f,men);n=n+1;
    VALID_L_1000_tobi = (VALID_L_1000_tobi>=1);
    VALID_R_1000_tobi = (VALID_R_1000_tobi>=1);

   
    
    %---------------------------------------------------------
%     
%       figure,plot(EYE_POS_LOG.pos_time,EYE_POS_LOG.pos_x), hold on 
%       figure,plot(time_x,POS_X_1000 )%  
%        plot(time_x(data_inter),POS_X_1000(data_inter),'o' )
%        plot(EYE_POS_LOG.pos_time)
%        
%        
%        figure
%         figure,plot(time_x_tobi,POS_X_1000_tobi) 
%        
    %---------------------------------------------------------
           
    % deletc extrem non-sense measure from TOBI
    extrem = abs(POS_X_1000_tobi)>3000;

    largo_log = numel(time_x);
    largo_tobi = numel(time_x_tobi);

    % synchronization 
    % search lag with minimun  error using the lag with max autocorrealtion! 

    % center data 
    C_POS_X_1000_tobi       = POS_X_1000_tobi - nanmedian(POS_X_1000_tobi);
    C_POS_X_1000            = POS_X_1000 - nanmedian(POS_X_1000);
    %C_POS_X_1000(data_inter)= 0;
    C_POS_X_1000(isnan(C_POS_X_1000))= 0;
    C_POS_X_1000_tobi(isnan(C_POS_X_1000_tobi)) = 0;
    C_POS_X_1000_tobi(extrem) = 0;

    C_POS_Y_1000_tobi       = POS_Y_1000_tobi - nanmedian(POS_Y_1000_tobi);
    C_POS_Y_1000            = POS_Y_1000 - nanmedian(POS_Y_1000);
    %C_POS_Y_1000(data_inter)= 0;
    C_POS_Y_1000(isnan(C_POS_Y_1000))= 0;
    C_POS_Y_1000_tobi(isnan(C_POS_Y_1000_tobi)) = 0;
    C_POS_Y_1000_tobi(extrem) = 0;

    % maximun evaluable Lag for autocorrelation 
    max_del = (largo_tobi-largo_log);
    
 
    
    % for the signal in X posicion 
    [Rxx, lagsx] = xcorr(C_POS_X_1000_tobi, C_POS_X_1000,max_del);%switch order of x and y
   % Rxx(lagsx<0)=0;
    [X, Ix] = max(Rxx);
    LAGx = lagsx(Ix);
    % for the signal in Y posicion 
    [Ryy, lagsy] = xcorr(C_POS_Y_1000_tobi, C_POS_Y_1000,max_del);%switch order of x and y
    %Ryy(lagsy<100000)=0;
    [Y, Iy] = min(Ryy); % NO SE PORQUE UNO QUEDO INVERTIDO ??
    LAGy = lagsy(Iy);

    [XY, Ixy] = max(-Ryy + Rxx); % NO SE PORQUE UNO QUEDO INVERTIDO ??
    LAGxy = lagsy(Ixy);

    disp([LAGx  LAGy LAGxy ] ) % set 69170 
    disp([X  -Y XY ] )% these values should be equal for (almost) perfect synchronization !

    % find the best synchronization with almos ALBITRARY criteria :) ...
    if sum([LAGx  LAGy LAGxy ]==mode([LAGx  LAGy LAGxy ])  )>1 % two equal lag .. this is good 
       LAG = mode([LAGx  LAGy LAGxy ])+1; 
    else   

    if std([LAGx  LAGy LAGxy ])>40; % low variability between lag .. the threshold no have any eviodence ...
        error('CHECK the synchronization !!!')
    else
        LAG = fix(mean([LAGx  LAGy LAGxy ]))+1;
    end

    end

%     %---------------------------------------------------------------------
%             % use  this line if there no well synchronization , to evaluate
%             % what happen 
               LAG = LAGx; %LAGy%   100;%
%             % % % CHECK the synchronization 
              figure,%_lan
              plot(C_POS_X_1000_tobi(LAG:LAG+largo_log)-nanmedian(C_POS_X_1000_tobi(LAG:LAG+largo_log))), hold on 
              plot(C_POS_X_1000); %
              figure,%_lan
              plot(-(C_POS_Y_1000_tobi(LAG:LAG+largo_log)-nanmedian(C_POS_Y_1000_tobi(LAG:LAG+largo_log)))), hold on 
              plot(C_POS_Y_1000)
              
              figure_lan
               plot(C_POS_X_1000_tobi)
                figure_lan
               plot(C_POS_Y_1000)
%     %---------------------------------------------------------------------
    
   end
    % convert no measure in NaN in .log data ... convertir dato no recogido a
    % NaN
    data_inter =  ~ismember(time_x,EYE_POS_LOG.pos_time);
    POS_X_1000(data_inter) = NaN;
    POS_Y_1000(data_inter) = NaN;
    data_inter_p =  ~ismember(time_pupil,EYE_PUPIL_LOG.pos_time);
    PUPIL_X_1000(data_inter) = NaN;
    PUPIL_Y_1000(data_inter) = NaN;

end
%% STEP:6
% add EYE tracker data to LAN structure 
if EYE
    syn_log  = time_x(1);
    if TOBI
    syn_tobi = time_x_tobi(LAG);

    presentation_time_tobi = time_x_tobi - syn_tobi + syn_log;
    end
    tt = timelan(LAN);

    nchan = 64;%LAN.nbchan;

    % add nan for the last poar of .loggf eye recording to abode matrix leght problem 
    additional_pionts=400000; 
    POS_X_1000(end+1:end+additional_pionts)=nan;
    POS_Y_1000(end+1:end+additional_pionts)=nan;
    PUPIL_X_1000(end+1:end+additional_pionts)=nan;
    PUPIL_Y_1000(end+1:end+additional_pionts)=nan;
    time_x(end+1:end+additional_pionts)=(time_x(end)+1):(time_x(end)+additional_pionts);


    for t = 1:length(LAN.data)
        % t=1
        if TOBI
        zero_tobi = find(presentation_time_tobi==fix(LAN.RT.OTHER.log_latency(t)));
        X_tobi = (zero_tobi+(LAN.time(t,1)*LAN.srate)):(zero_tobi+(LAN.time(t,2)*LAN.srate));
        X_tobi =  X_tobi(1:end-1);
        end
        zero_log = find(time_x==fix(LAN.RT.OTHER.log_latency(t)));% ajustando los tiempos del EEG con los del eye-tracker
        X_log = (zero_log+(LAN.time(t,1)*LAN.srate)):(zero_log+(LAN.time(t,2)*LAN.srate));
        X_log =  X_log(1:end-1);

        red=3;
        amp=30;
        if TOBI
        LAN.data{t}(nchan+1,:) = POS_X_1000_tobi(X_tobi)/red;
        
        LAN.data{t}(nchan+3,:) = POS_Y_1000_tobi(X_tobi)/red;
        

        LAN.data{t}(nchan+5,:) = POS_XL_1000_tobi(X_tobi)/red;
        LAN.data{t}(nchan+6,:) = POS_YL_1000_tobi(X_tobi)/red;
        LAN.data{t}(nchan+7,:) = POS_XR_1000_tobi(X_tobi)/red;
        LAN.data{t}(nchan+8,:) = POS_YR_1000_tobi(X_tobi)/red;

        
       

        LAN.data{t}(nchan+9,:) = PUPIL_R_1000_tobi(X_tobi)*amp;
        
        LAN.data{t}(nchan+11,:) = PUPIL_L_1000_tobi(X_tobi)*amp;
        

        LAN.data{t}(nchan+13,:) = Gdir_XL_1000_tobi(X_tobi);
        LAN.data{t}(nchan+14,:) = Gdir_YL_1000_tobi(X_tobi);
        LAN.data{t}(nchan+15,:) = Gdir_ZL_1000_tobi(X_tobi);

        LAN.data{t}(nchan+16,:) = Gdir_XR_1000_tobi(X_tobi);
        LAN.data{t}(nchan+17,:) = Gdir_YR_1000_tobi(X_tobi);
        LAN.data{t}(nchan+18,:) = Gdir_ZR_1000_tobi(X_tobi);

        LAN.data{t}(nchan+19,:) = VALID_L_1000_tobi(X_tobi); 
        LAN.data{t}(nchan+20,:) = VALID_R_1000_tobi(X_tobi); 
        
        end
         
        LAN.data{t}(nchan+2,:) = POS_X_1000(X_log)/red;
        LAN.data{t}(nchan+4,:) = POS_Y_1000(X_log)/red;
        LAN.data{t}(nchan+10,:) = PUPIL_Y_1000(X_log)*amp;
        LAN.data{t}(nchan+12,:) = PUPIL_X_1000(X_log)*amp;
    end


    LAN.chanlocs(nchan+1).labels = 'X_tobi' ; LAN.chanlocs(nchan+1).type = 'EYE' ;
    LAN.chanlocs(nchan+2).labels = 'X_log' ; LAN.chanlocs(nchan+2).type = 'EYE' ;
    LAN.chanlocs(nchan+3).labels = 'Y_tobi' ; LAN.chanlocs(nchan+3).type = 'EYE' ;
    LAN.chanlocs(nchan+4).labels = 'Y_log' ; LAN.chanlocs(nchan+4).type = 'EYE' ;

    LAN.chanlocs(nchan+5).labels = 'XL_tobi' ; LAN.chanlocs(nchan+5).type = 'EYE' ;
    LAN.chanlocs(nchan+6).labels = 'YL_log' ; LAN.chanlocs(nchan+6).type = 'EYE' ;
    LAN.chanlocs(nchan+7).labels = 'XR_tobi' ; LAN.chanlocs(nchan+7).type = 'EYE' ;
    LAN.chanlocs(nchan+8).labels = 'YR_log' ; LAN.chanlocs(nchan+8).type = 'EYE' ;

    LAN.chanlocs(nchan+9).labels = 'P_R_tobi' ; LAN.chanlocs(nchan+9).type = 'EYE' ;
    LAN.chanlocs(nchan+10).labels = 'P_R_log' ; LAN.chanlocs(nchan+10).type = 'EYE' ; % Y is Rigth 
    LAN.chanlocs(nchan+11).labels = 'P_L_tobi' ; LAN.chanlocs(nchan+11).type = 'EYE' ;
    LAN.chanlocs(nchan+12).labels = 'P_L_log' ; LAN.chanlocs(nchan+12).type = 'EYE' ; % X is left 

    LAN.chanlocs(nchan+13).labels = 'GD_XL_tobi' ; LAN.chanlocs(nchan+13).type = 'EYE' ;
    LAN.chanlocs(nchan+14).labels = 'GD_YL_tobi' ; LAN.chanlocs(nchan+14).type = 'EYE' ;
    LAN.chanlocs(nchan+15).labels = 'GD_ZL_tobi' ; LAN.chanlocs(nchan+15).type = 'EYE' ;

    LAN.chanlocs(nchan+16).labels = 'GD_XR_tobi' ; LAN.chanlocs(nchan+16).type = 'EYE' ;
    LAN.chanlocs(nchan+17).labels = 'GD_YR_tobi' ; LAN.chanlocs(nchan+17).type = 'EYE' ;
    LAN.chanlocs(nchan+18).labels = 'GD_ZR_tobi' ; LAN.chanlocs(nchan+18).type = 'EYE' ;

    LAN.chanlocs(nchan+19).labels = 'Valid_L_tobi' ; LAN.chanlocs(nchan+19).type = 'EYE' ;
    LAN.chanlocs(nchan+20).labels = 'Valid_R_tobi' ; LAN.chanlocs(nchan+20).type = 'EYE' ;


    LAN = lan_check(LAN);

    %prepro_plot(LAN);

  

end

  %save current version of this script in the LAN strucutre (guarda el script del procesamiento del dato específico)

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
  
% % CHECK synchronization : SEE: paterns of blink and missing data 
% %                              and... congruence between  log and tobii data    
% prepro_plot(LAN)

%%
% save in local disck
%save([ PATH_D 'MAT/LAN_' TAREA '_array_EYE'],'LAN', '-v7.3');
mkdir([  DATA '_pro' ])
save([DATA '_pro/LAN_' TAREA '_array_EYE' ] , 'LAN', '-v7.3')


%save in neuroCICS DB

% mkdir(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG' ] , [  DATA '_pro' ])
% save(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG/' DATA '_pro/LAN_' TAREA '_array_EYE' ] , 'LAN', '-v7.3')

%clear

toc % 
