% 
% Model RL from 
% Time-frequiency amplitude 
%


clear
% eeg analisis 
Sujetos = {{
%'CFOI_23101987' % 'RL WM GN'     %
%'CAMS_24101973' % 'RL WM GN'  %
%'FLFS_08081990' % 'RL GN WM'  % 
%'GAMY_03071978' % RL GN WM'  % 
%'K_BA_22111991' % RL GN WM % 
%'KSDG_20051983' % RL GN WM %  
'MMPM_06011977' % 'RL GN WM' 
%'S_WL_25051979' % 'RL GN WM'
%'CAPG_27061969' % 'RL GN WM'
%'MJDD_01101984'
}',{% 'RL GN WM' 
'CEAS_23071992' % 'RL GN WM' 
%'CAHG_27061988' % 'RL GN WM' 
%'YDCL_11021994' % 'RL GN WM' 
%'B_FE_24011989' % 'RL GN WM' 
%'MEFR_14061991' % 'RL GN WM' 
%'SNCS_22121989' % 'RL GN WM' 
%'MAAO_15081988' % 'RL GN WM'  
%'GRCN_21111979' % 'RL GN WM'  
%'KXMS_13061971' % 'RL GN WM' 
%'CFPL_08021984' % 'RL GN WM' 
}'} ;

GRS={'pacientes','controles'};

OVERWRITE=false;



for nG = 1:2
    GR= GRS{nG};
    Sujetos_g = Sujetos{nG};
for  nS =Sujetos_g %%
    
for nT = {'WM'}
   disp(['Sujetos ' nS{1} ' task: ' nT{1}  ])
   
   
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
    elseif ~OVERWRITE &&  exist ([ PATH_MAT 'LAN_' TAREA '_MODEL_ML_SMP.mat' ], 'file')
      fprintf(['OK algunos modelos  ' SU ' se revisara cada uno .... ojo:)  \n'])  
     % continue
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
  

   % Correct responses 
   % Succefull Memory performance 
   nB1=90;
   nB2=90;
   if LAN.trials<180 && LAN.trials>90
       %error('fixme')
       %continue
       nB1=LAN.trials-90;
   elseif LAN.trials==90
       nB1=90;
       nB2=0;
   elseif LAN.trials<90
       %error('fixme')
       %continue 
       nB1=LAN.trials;
       nB2=0;
   end
   % check responces per  bloq 
   B1 =  [ ones([1 nB1])  zeros([ 1 nB2])  ];
   b1r20 = LAN.RT.resp(LAN.RT.OTHER.log_probe==20 & B1 );
   b1r21 = LAN.RT.resp(LAN.RT.OTHER.log_probe==21 & B1 );
   b1r40 = LAN.RT.resp(LAN.RT.OTHER.log_probe==40 & B1 );
   b1r41 = LAN.RT.resp(LAN.RT.OTHER.log_probe==41 & B1 );
   b1r60 = LAN.RT.resp(LAN.RT.OTHER.log_probe==60 & B1 );
   b1r61 = LAN.RT.resp(LAN.RT.OTHER.log_probe==61 & B1 );
   if mean( [ b1r20==1  b1r21==2  ]) > 0.5 || mean( [ b1r40==1  b1r41==2  ])>0.5
   B1c=2;
   elseif mean( [ b1r20==2  b1r21==1  ]) > 0.5 || mean( [ b1r40==2  b1r41==1  ])>0.5
   B1c=1;    
   elseif mean( [  b1r21==2   b1r41==2 ]) > 0.5 
   B1c=2;
   elseif mean( [  b1r21==1   b1r41==1 ]) > 0.5 
   B1c=1;    
   else
       B1c=1;  
       %error('Respuestas no consistentes')
   end
   
   B2 = ~B1;
   b2r20 = LAN.RT.resp(LAN.RT.OTHER.log_probe==20 & B2 );
   b2r21 = LAN.RT.resp(LAN.RT.OTHER.log_probe==21 & B2);
   b2r40 = LAN.RT.resp(LAN.RT.OTHER.log_probe==40 & B2 );
   b2r41 = LAN.RT.resp(LAN.RT.OTHER.log_probe==41 & B2);
   b2r60 = LAN.RT.resp(LAN.RT.OTHER.log_probe==60 & B2 );
   b2r61 = LAN.RT.resp(LAN.RT.OTHER.log_probe==61 & B2);
   if mean( [ b2r20==1  b2r21==2  ]) > 0.5 || mean( [ b2r40==1  b2r41==2  ])>0.5
      % 
   B2c=2;
   elseif mean( [ b2r20==2  b2r21==1 b2r40==2  b2r41==1   ]) > 0.5 || mean( [ b2r40==2  b2r41==1   ])>0.5
   B2c=1;   
   elseif mean( [  b2r21==2   b2r41==2 ]) > 0.5 
   B2c=2;
   elseif mean( [  b2r21==1   b2r41==1 ]) > 0.5 
   B2c=1;    
   elseif B2c==1;  
          B2c=2;
   elseif B2c==2;  
          B2c=1;          
   else       
       error('Respuestas no consistentes')
   end 
   
   if B2c == B1c
       warning('Parece qie sujeto no cambio de mano')
   end
   
   if B1c==1 
      B1c0=2;
   else
      B1c0=1;
   end
   if B2c==1 
      B2c0=2;
   else
      B2c0=1;
   end
   
   
   fprintf(['-------PERFORMNACE--------\n']);
   fprintf(['\t\t B1 \t\t B2\n']);
   fprintf(['load2\t\t '   num2str(mean( [b1r20==B1c0    b1r21==B1c ]),3)            '\t\t'    num2str(mean( [b2r20==B2c0    b2r21==B2c ]),3)        '\n']);
   fprintf(['load4\t\t '   num2str(mean( [b1r40==B1c0    b1r41==B1c ]),3)            '\t\t'    num2str(mean( [b2r40==B2c0    b2r41==B2c ]),3)        '\n']);
   fprintf(['load6\t\t '   num2str(mean( [b1r60==B1c0    b1r61==B1c ]),3)            '\t\t'    num2str(mean( [b2r60==B2c0    b2r61==B2c ]),3)        '\n']);
   fprintf(['------------------------\n']);
   
   
    SMP_b1_1 = (LAN.RT.OTHER.log_probe==21|LAN.RT.OTHER.log_probe==41|LAN.RT.OTHER.log_probe==61) & (B1c==LAN.RT.resp); 
    SMP_b1_0 = (LAN.RT.OTHER.log_probe==20|LAN.RT.OTHER.log_probe==40|LAN.RT.OTHER.log_probe==60) & (B1c0==LAN.RT.resp); 
    SMP_b1 = SMP_b1_1 +SMP_b1_0;
    
    SMP_b2_1 = (LAN.RT.OTHER.log_probe==21|LAN.RT.OTHER.log_probe==41|LAN.RT.OTHER.log_probe==61) & (B1c==LAN.RT.resp); 
    SMP_b2_0 = (LAN.RT.OTHER.log_probe==20|LAN.RT.OTHER.log_probe==40|LAN.RT.OTHER.log_probe==60) & (B1c0==LAN.RT.resp); 
    SMP_b2 = SMP_b2_1 + SMP_b2_0;
    
    SMP  = [SMP_b1(1:nB1)  SMP_b2(nB1+1:end)];
    LOAD = (LAN.RT.OTHER.log_code==12)*2 + (LAN.RT.OTHER.log_code==14)*4 + (LAN.RT.OTHER.log_code==16)*6 ;
    % time
    T1 = find_approx(LAN.freq.time, -0.5);
    T2 = find_approx(LAN.freq.time, 3.8); % encoiding and mantenances 
   
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
         if ~OVERWRITE &&  exist ([ PATH_MAT 'LAN_' TAREA '_MODEL_ML_SMP.mat' ], 'file')
             disp([SU TAREA ' M1 OK'])
         else
         cfgM=[];
         cfgM.type = 'glm';
         cfgM.ops = ['pre(glm:_M1_' SU '_' TAREA  '  )'];
         [w, ww] =  lan_model_stat(FT(:,:,T1:T2,:),ones(size(LOAD(~noLAN))),LOAD(~noLAN), SMP(~noLAN),cfgM);
         
         LAN.freq.model.t =ww.t;
         LAN.freq.model.p= w;
         LAN.freq.model.b= ww.b;
         LAN.freq.model.r= {'Int' ,'LOAD','SMP'};
         LAN.freq.model.freq=LAN.freq.freq;
         LAN.freq.model.time=LAN.freq.time(T1:T2);
         
         LAN.data=[];
         
         save ([ PATH_MAT 'LAN_' TAREA '_MODEL_ML_SMP' ],'LAN', '-v7.3')
         ntry=0;
         while ntry<5
             try         
         %save(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG/' DATE '_pro/LAN_' TAREA '_MODEL_ML_SMP' ] , 'LAN', '-v7.3')
                ntry=10;
             catch
                pause(rand(1)*30)
                ntry=ntry+1;
             end
         end  
         end
         
         
         if ~OVERWRITE &&  exist ([ PATH_MAT 'LAN_' TAREA '_MODEL_ML_SMP_noise2.mat' ], 'file')
             disp([SU TAREA '  M2 OK'])
         else
         cfgM=[];
         cfgM.type = 'glm';
         cfgM.ops = ['pre(glm:_M2_' SU '_' TAREA ' + Noise )'];
         [w, ww] =  lan_model_stat(FT(:,:,T1:T2,:),ones(size(LOAD(~noLAN))), LOAD(~noLAN), SMP(~noLAN),NoiseLevel(~noLAN),LowSignal(~noLAN),cfgM);         
         LAN.freq.model.t =ww.t;
         LAN.freq.model.p= w;
         LAN.freq.model.b= ww.b;
         LAN.freq.model.r= {'Int' ,'LOAD','SMP','noiseNL', 'noiseLS'};
         LAN.freq.model.freq=LAN.freq.freq;
         LAN.freq.model.time=LAN.freq.time(T1:T2);
         
         LAN.data=[];
    
         save ([ PATH_MAT 'LAN_' TAREA '_MODEL_ML_SMP_noise2' ],'LAN', '-v7.3')
         ntry=0;
         while ntry<5
             try         
         %save(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG/' DATE '_pro/LAN_' TAREA '_MODEL_ML_SMP_noise2' ] , 'LAN', '-v7.3')
                ntry=10;
             catch
                pause(rand(1)*30)
                ntry=ntry+1;
             end
         end  
         end
         
         
         
         if ~OVERWRITE &&  exist ([ PATH_MAT 'LAN_' TAREA '_MODEL_ML_SMP_MLxSMP.mat' ], 'file')
             disp([SU TAREA '  M3 OK'])
         else
         cfgM=[];
         cfgM.type = 'glm';
         cfgM.ops = ['pre(glm:_M3_' SU '_' TAREA '  )'];
         [w, ww] =  lan_model_stat(FT(:,:,T1:T2,:),ones(size(LOAD(~noLAN))),LOAD(~noLAN), SMP(~noLAN), LOAD(~noLAN).*SMP(~noLAN), cfgM);
         
         LAN.freq.model.t =ww.t;
         LAN.freq.model.p= w;
         LAN.freq.model.b= ww.b;
         LAN.freq.model.r= {'Int' ,'ML','SMP','ML*SMP'};
         LAN.freq.model.freq=LAN.freq.freq;
         LAN.freq.model.time=LAN.freq.time(T1:T2);
         
         LAN.data=[];
         
         save ([ PATH_MAT 'LAN_' TAREA '_MODEL_ML_SMP_MLxSMP' ],'LAN', '-v7.3')
         ntry=0;
         while ntry<5
             try         
         %save(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG/' DATE '_pro/LAN_' TAREA '_MODEL_ML_SMP_MLxSMP' ] , 'LAN', '-v7.3')
                ntry=10;
             catch
                pause(rand(1)*30)
                ntry=ntry+1;
             end
         end
         end
         
         
         if ~OVERWRITE &&  exist ([ PATH_MAT 'LAN_' TAREA '_MODEL_ML_SMP_MLxSMP_noise2.mat' ], 'file')
             disp([SU TAREA 'M4 OK'])
         else
         cfgM=[];
         cfgM.type = 'glm';
         cfgM.ops = ['pre(glm:_M4_' SU '_' TAREA  ' + Noise )'];
         [w, ww] =  lan_model_stat(FT(:,:,T1:T2,:),ones(size(LOAD(~noLAN))),LOAD(~noLAN), SMP(~noLAN), LOAD(~noLAN).*SMP(~noLAN),NoiseLevel(~noLAN),LowSignal(~noLAN),cfgM);         
         LAN.freq.model.t =ww.t;
         LAN.freq.model.p= w;
         LAN.freq.model.b= ww.b;
         LAN.freq.model.r= {'Int' ,'ML','SMP','ML*SMP','noiseNL', 'noiseLS'};
         LAN.freq.model.freq=LAN.freq.freq;
         LAN.freq.model.time=LAN.freq.time(T1:T2);
         
         LAN.data=[];
    
         save ([ PATH_MAT 'LAN_' TAREA '_MODEL_ML_SMP_MLxSMP_noise2' ],'LAN', '-v7.3')
         ntry=0;
         while ntry<5
             try
                %save(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG/' DATE '_pro/LAN_' TAREA '_MODEL_ML_SMP_MLxSMP_noise2' ] , 'LAN', '-v7.3')
                ntry=10;
             catch
                pause(rand(1)*30)
                ntry=ntry+1;
             end
         end
         end
    






end
end
end


% %% for checking data de data is OK
if 0 
figure
pcolor2(LAN.freq.model.time, LAN.freq.model.freq, squeeze(LAN.freq.model.t{4}(:,59,:))), shading flat
end


