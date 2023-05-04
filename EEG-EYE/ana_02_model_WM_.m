% 
% Model RL from 
% Time-frequiency amplitude 
%


clear
% eeg analisis 
Sujetos = {
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
'GADM_08111983' % 'RL GN WM' % sp1 
'GAQJ_21061987' % 'RL GN WM' 
'IABF_16091997' % 'RL GN WM' 
'IABM_03061982' % 'RL GN WM' %RRR
'JAPV_03111995' % 'RL GN WM'
'JAZV_27081991' % 'RL GN WM'
'JCZD_03071963' % 'RL GN WM'% RRR % sp1 
'JDRP_06081956' % 'RL GN WM'  % sp1 
'JFFA_04101976' % 'RL GN WM'  % sp1 
'JFSO_05061967' % 'RL GN WM'  % sp1 
'JGRF_24091991' % 'RL GN WM'  % sp1 
 'LAAB_29101982' % 'RL GN WM'  % sp1
'LDEV_31011988' % 'RL GN WM'  % sp1
'LMRG_21041962' % 'RL GN WM'  % sp1
'MAEM_17051957' % 'RL GN WM'  % sp1
'MFCM_21011984' % 'RL GN WM'  % sp1
'MGRN_22101979' % 'RL GN WM'  % sp1
'MLAD_09021957' % 'RL GN WM'  % sp1
'MLCB_11041981' % 'RL GN WM'  % sp1
'MPMM_22111977' % 'RL GN WM'  % sp1
'MTVV_27111974' % 'RL GN WM'  % sp1
'OAPC_01041962' % 'RL GN WM'  % sp1
'RABG_09091981' % 'RL GN WM' 
'RAPA_19121977' % 'RL GN WM'  % sp1
'REOO_10031989' % 'RL GN WM'  % sp1
'RFCR_12111978' % 'RL GN WM'  % sp1
'RSVT_17091955' % 'RL GN WM'  % sp1
'SAAA_10051991' % 'RL GN WM'  % sp1
'SAJC_22021972' % 'RL GN WM'  % sp1
'SAVS_23061979' % 'RL GN WM'  % sp1
'UENM_22121979' % 'RL GN WM'  % sp1
'VAPH_11101985' % 'RL GN WM'  % sp1
}' ;


OVERWRITE=false;


for  nS =Sujetos(1:end) %%
for nT = {'WM'}
  
    SU=     nS{1}; 
    disp(SU)
    PATH_MAT =  '/Volumes/DBNC_03/neuroCOVID/DATA/%S/EEG/%D_pro/';
    PATH_EEG =  '/Volumes/DBNC_03/neuroCOVID/DATA/%S/EEG/';
    PATH_LOG =  '/Volumes/DBNC_03/neuroCOVID/DATA/%S/LOG/';
    
     
    GR = '';
    TAREA = nT{1};       %
    DATE = '';


    PATH_EEG =strrep( strrep(PATH_EEG,'%G',GR) , '%S' , SU);
    if isempty(DATE)
       [paso paso2] = unix(['ls  -d ' PATH_EEG '*_pro' ]);
        DATE = paso2(end-12:end-5);
    end
    PATH_MAT = strrep( strrep( strrep(PATH_MAT,'%G',GR) , '%S' , SU), '%D', DATE);
    PATH_LOG = strrep( strrep( strrep(PATH_LOG,'%G',GR) , '%S' , SU), '%D', DATE);
    
    
    if ~exist ([ PATH_MAT 'LAN_' TAREA '_EEG_interp_freq.mat' ], 'file')
      fprintf(['No data  ' SU ' \n'])  
      continue
    elseif ~OVERWRITE &&  exist ([ PATH_MAT 'LAN_' TAREA '_MODEL_ML_SMP.mat' ], 'file')
      fprintf(['OK algunos modelos  ' SU ' \n'])  
      %continue
    end

    cd([PATH_MAT ])
    load ([ PATH_MAT 'LAN_' TAREA '_EEG_interp_freq' ])
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
         save(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG/' DATE '_pro/LAN_' TAREA '_MODEL_ML_SMP' ] , 'LAN', '-v7.3')
                ntry=10;
             catch
                pause(rand(1)*30)
                ntry=ntry+1;
             end
         end  
         end
         
         
         if ~OVERWRITE &&  exist ([ PATH_MAT 'LAN_' TAREA '_MODEL_ML_SMP_noise2.mat' ], 'file')
             disp([SU TAREA '  M1 OK'])
         else
         cfgM=[];
         cfgM.type = 'glm';
         cfgM.ops = ['pre(glm:_M1_' SU '_' TAREA ' + Noise )'];
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
         save(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG/' DATE '_pro/LAN_' TAREA '_MODEL_ML_SMP_noise2' ] , 'LAN', '-v7.3')
                ntry=10;
             catch
                pause(rand(1)*30)
                ntry=ntry+1;
             end
         end  
         end
         
         
         
         if ~OVERWRITE &&  exist ([ PATH_MAT 'LAN_' TAREA '_MODEL_ML_SMP_MLxSMP.mat' ], 'file')
             disp([SU TAREA '  M2 OK'])
         else
         cfgM=[];
         cfgM.type = 'glm';
         cfgM.ops = ['pre(glm:_M2_' SU '_' TAREA '  )'];
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
         save(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG/' DATE '_pro/LAN_' TAREA '_MODEL_ML_SMP_MLxSMP' ] , 'LAN', '-v7.3')
                ntry=10;
             catch
                pause(rand(1)*30)
                ntry=ntry+1;
             end
         end
         end
         
         
         if ~OVERWRITE &&  exist ([ PATH_MAT 'LAN_' TAREA '_MODEL_ML_SMP_MLxSMP_noise2.mat' ], 'file')
             disp([SU TAREA 'M2 OK'])
         else
         cfgM=[];
         cfgM.type = 'glm';
         cfgM.ops = ['pre(glm:_M2_' SU '_' TAREA  ' + Noise )'];
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
                save(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG/' DATE '_pro/LAN_' TAREA '_MODEL_ML_SMP_MLxSMP_noise2' ] , 'LAN', '-v7.3')
                ntry=10;
             catch
                pause(rand(1)*30)
                ntry=ntry+1;
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


