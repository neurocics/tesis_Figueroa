% 
% Model RL from 
% Time-frequiency amplitude 
%

% Time-frequiency calculations 
clear
% eeg analisis 
Sujetos = {{
%'CFOI_23101987' % 'RL WM GN'     %
%'CAMS_24101973' % 'RL WM GN'  %
%'FLFS_08081990' % 'RL GN WM'  % 
%'GAMY_03071978' % RL GN WM'  % 
%'K_BA_22111991' % RL GN WM % 
%'KSDG_20051983' % RL GN WM %  
%'MMPM_06011977' % 'RL GN WM' 
'S_WL_25051979' % 'RL GN WM'
'CAPG_27061969' % 'RL GN WM'
%'MJDD_01101984'
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

OVERWRITE=true;%false;

PATH_MAT =  '/Volumes/Alehermosa/TESIS/%S/EEG/%D_pro/';
PATH_EEG =  '/Volumes/Alehermosa/TESIS/%S/EEG/';




for nG = 1:2
    GR= GRS{nG};
    Sujetos_g = Sujetos{nG};
for  nS =Sujetos_g %%
    
for nT = {'GN'}
    disp(['Sujeto ' nS{1} ' task: ' nT{1}  ])
   
   
    PATH_MAT =  '/Volumes/Alehermosa/TESIS/%G/%S/EEG/%D_pro/';
    PATH_EEG =  '/Volumes/Alehermosa/TESIS/%G/%S/EEG/';
    PATH_LOG =  '/Volumes/Alehermosa/TESIS/%G/%S/LOG/';
    
    SU=     nS{1};  
    TAREA = nT{1};       %
    DATE = '';


    PATH_EEG =strrep( strrep(PATH_EEG,'%G',GR) , '%S' , SU);
    if isempty(DATE)
       [paso paso2] = unix(['ls  -d ' PATH_EEG '*_pro' ]);
        DATE = paso2(end-12:end-5);
    end
    PATH_MAT = strrep( strrep( strrep(PATH_MAT,'%G',GR) , '%S' , SU), '%D', DATE);
    PATH_LOG = strrep( strrep( strrep(PATH_LOG,'%G',GR) , '%S' , SU), '%D', DATE);
    
    
    if ~exist ([ PATH_MAT 'LAN_' TAREA '_EEG_interp_freq_lapla.mat' ], 'file')
      fprintf(['No data  ' SU ' \n'])  
      continue
    elseif ~OVERWRITE &&  exist ([ PATH_MAT 'LAN_' TAREA '_MODEL_NGO_nGO_nNGP_G_pG_Lf_noise2.mat' ], 'file')
      fprintf(['modelo listos  OK ' SU ' \n'])  
      continue
    end

    cd([PATH_MAT ])
    load ([ PATH_MAT 'LAN_' TAREA '_EEG_interp_freq_lapla' ])
    FT = lan_getdatafile(LAN.freq.powspctrm.filename,LAN.freq.powspctrm.path,LAN.freq.powspctrm.trials);


    noFT= isemptycell(FT);
    noLAN = LAN.accept==0;
    if any(noFT~=noLAN)
        error('Data no calsa!')
    else
        FT =cat(4,FT{:});
    end
  
    if isfield(LAN.RT, 'OTHER') && isfield(LAN.RT.OTHER, 'log_code')
        GO = LAN.RT.OTHER.log_code==10;
        NGO = LAN.RT.OTHER.log_code==21;
    else
        GO  = LAN.RT.est==4 | LAN.RT.est==10;
        NGO = LAN.RT.est==5 | LAN.RT.est==21;
    end
    
    % Check the stimuli
    if sum(NGO)*2 >= sum(GO)
           error('estimulos incosistentes')
    end
    
    % sequences 
    nGO=ones(1,LAN.trials);
    for t = 2:LAN.trials
    if GO(t)==1
        nGO(t) = nGO(t-1)+1;
    elseif NGO(t)==1
        nGO(t) = 0;   
    else
        error('estimulos incosistentes')
    end
    end
    nNGO = zeros(1,LAN.trials);
    nNGO(nGO==0) = nGO(find(nGO==0)-1);
    % check this! % trucate effect 
    %nNGO(nNGO>5)=6;
    %nGO(nGO>5)=6;
   % check responces per  bloq 
Q=0.25;
seq=[];
seq(1)=1;
for t = 2:length(LAN.RT.est)
   if LAN.RT.est(t)~=LAN.RT.est(t-1)
       seq(t)=1;
   else
       seq(t)=seq(t-1)+1;
   end

end

Seq=seq;
%UNO=Seq==1 & LAN.RT.est==2;
ind= find(LAN.RT.est==5);
ind=ind(ind>1);
%Seq(ind) = Seq(ind-1)+1;
Seq = 1- ((1-Q).^(Seq-1));
Seq(ind) = 1-Seq(ind-1) ;
   
   %correct
   if sum(LAN.RT.resp>0)==0
      %try 
           R = ls_lan([PATH_LOG  SU '_' DATE '_GN-gonogo_*lum*.log']);
           ind = isnan(fun_in_cell(R,'(strfind(@,''training''))')) ; 
           if sum(ind)~=1
               disp(['Revisar GN >>>>>>>>  ' SU  ])
               %continue %% ver los saltados 

           end
           R = R{ind} ;
                cfg             = [];
                cfg.filename    = R;
                cfg.type        = 'presentation';
                cfg.est         = [10 21];
                cfg.resp        = [1 2 3 4 5 6] ;
                cfg.rw          = [5000] ;      %   (ms)    %   (ms)
                RTlog = rt_read(cfg);

                if length(LAN.RT.est) ~= length(RTlog.est)
                    error('numero respuestas  incosistentes')
                end
      %catch
           % lan_read_file()        
          
      %end
        
   end
   
   GOOD_GO = (RTlog.est==10)&(RTlog.resp>0);
   GOOD_NGO = (RTlog.est==21)&(RTlog.resp==-99);
   GOOD =  GOOD_GO + GOOD_NGO;
   pGOOD = [ 0 GOOD(1:end-1)];
   
   LEFT=RTlog.resp==1;
   
   fprintf(['-------PERFORMNACE--------\n']);
   fprintf(['\t\t B1 \t\n']);
   fprintf(['GO\t\t '    num2str(mean(GOOD_GO(GO)),3)       '\n']);
   fprintf(['NGO\t\t '   num2str(mean(GOOD_NGO(NGO)),3)    '\n']);
   fprintf(['------------------------\n']);
   
   
     % time
    T1 = find_approx(LAN.freq.time, -0.5);
    T2 = find_approx(LAN.freq.time, 1); % encoiding and mantenances 
   
    %%% calcular regresores para el modelo 
    
    % Noise Nuance Regressor     
    if any (LAN.tag.mat(:)>0)
      LowSignal = mean(LAN.tag.mat>0,1) ;
    else
      LowSignal = load ([ PATH_MAT 'LAN_' TAREA '_EEG_interp' ]); 
      LowSignal = normal_z(mean(LowSignal.LAN.tag.mat>0,1)) ;
    end
    % Global high frequiency 
    H1 =find_approx(LAN.freq.fourierp.freq,15);
    H2 =find_approx(LAN.freq.fourierp.freq,45);
    NoiseLevel = normal_z( squeeze(nanmean(nanmean(LAN.freq.fourierp.data(:,H1:H2,:)        ,1),2))');
    
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
         

    
    %%% fitting model 
    
         cfgM=[];
         cfgM.type = 'glm';
         cfgM.ops = ['pre(glm:_M1_' SU '_' TAREA  '  )'];
         [w, ww] =  lan_model_stat(FT(:,:,T1:T2,:),ones(size(NGO(~noLAN))),....
             NGO(~noLAN), ...
             GO(~noLAN).*Seq(~noLAN)  ,NGO(~noLAN).*Seq(~noLAN), ...
             GOOD(~noLAN),pGOOD(~noLAN),....
             LEFT(~noLAN), cfgM);
         
         LAN.freq.model.t =ww.t;
         LAN.freq.model.p= w;
         LAN.freq.model.b= ww.b;
         LAN.freq.model.r= {'Int' ,'NGO','nGO','nNGO','GOOD','pGOOD','LEFT'};
         LAN.freq.model.freq=LAN.freq.freq;
         LAN.freq.model.time=LAN.freq.time(T1:T2);
         
         LAN.data=[];
         
         save ([ PATH_MAT 'LAN_' TAREA '_MODEL_NGO_nGO_nNGP_G_pG_Lf' ],'LAN', '-v7.3')
        % save(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG/' DATE '_pro/LAN_' TAREA '_MODEL_NGO_nGO_nNGP_G_pG_Lf'] , 'LAN', '-v7.3')

         cfgM=[];
         cfgM.type = 'glm';
         cfgM.ops = ['pre(glm:_M1_' SU '_' TAREA ' + Noise )'];
         [w, ww] =  lan_model_stat(FT(:,:,T1:T2,:),ones(size(NGO(~noLAN))),...
             NGO(~noLAN), ...nGO(~noLAN),...
             GO(~noLAN).*Seq(~noLAN), ....
             NGO(~noLAN).*Seq(~noLAN), ...
             GOOD(~noLAN),pGOOD(~noLAN),LEFT(~noLAN),NoiseLevel(~noLAN),LowSignal(~noLAN),cfgM);         
         LAN.freq.model.t =ww.t;
         LAN.freq.model.p= w;
         LAN.freq.model.b= ww.b;
         LAN.freq.model.r= {'Int' ,'NGO','nGO','nNGO','GOOD','pGOOD','LEFT','noiseNL', 'noiseLS'};
         LAN.freq.model.freq=LAN.freq.freq;
         LAN.freq.model.time=LAN.freq.time(T1:T2);
         
         LAN.data=[];
    
         save ([ PATH_MAT 'LAN_' TAREA '_MODEL_NGO_nGO_nNGP_G_pG_Lf_noise2' ],'LAN', '-v7.3')
         %save(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG/' DATE '_pro/LAN_' TAREA '_MODEL_NGO_nGO_nNGP_G_pG_Lf_noise2' ] , 'LAN', '-v7.3')
        
         

end
end
end

% %% for checking data de data is OK
if 0 
figure
pcolor2(LAN.freq.model.time, LAN.freq.model.freq, squeeze(LAN.freq.model.t{3}(:,65,:))), shading flat
end


