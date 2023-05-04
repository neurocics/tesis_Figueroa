% 
% Model RL from 
% Time-frequiency amplitude 
%


clear
% eeg analisis 
Sujetos = {
'AGBB_26121972' % 'RL WM GN'     % OK
'ASCS_31121976' % 'RL WM GN'  % sp1
'BBMA_19041977' % 'RL GN WM'  % sp1
'BCUV_02091982' % RL GN WM'  % sp1
'C_DG_15111995' % RL GN WM % sp1 GN
'CACB_27091965' % RL GN WM % sp1 
'CAMO_08111988' % 'RL GN WM' 
'CARM_14011974' % 'RL GN WM'
'CDPR_25031986' % 'RL GN WM'
'CEBT_30091989' % 'RL GN WM' % sp1 
'CJGV_24021957' % 'RL GN WM' % sp1 
'CMCO_17051975' % 'RL GN WM' 
'COGC_18022000' % 'RL GN WM' 
'EERV_01041955' % 'RL GN WM' 
'ELFN_09021961' % 'RL GN WM' 
'FABN_19111984' % 'RL GN WM' 
'FJAM_10012001' % 'RL GN WM' % sp1 
'FRCA_12081956' % 'RL GN WM' % sp1 
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
'UENM_22121979' % 'RL GN WM'  % sp1 falta LOG
'VAPH_11101985' % 'RL GN WM'  % sp1
}' ;


OVERWRITE=false;


for  nS =Sujetos %%
for nT = {'GN'}
   disp('Sujetos(1)')
   
   
    PATH_MAT =  '/Volumes/DBNC_03/neuroCOVID/DATA/%S/EEG/%D_pro/';
    PATH_EEG =  '/Volumes/DBNC_03/neuroCOVID/DATA/%S/EEG/';
    PATH_LOG =  '/Volumes/DBNC_03/neuroCOVID/DATA/%S/LOG/';
    
    SU=     nS{1};  
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
    elseif ~OVERWRITE &&  exist ([ PATH_MAT 'LAN_' TAREA '_MODEL_NGO_nGO_nNGP_G_pG_Lf_noise2.mat' ], 'file')
      fprintf(['freq OK ' SU ' \n'])  
      continue
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
    nNGO(nNGO>5)=6;
    nGO(nGO>5)=6;
   % check responces per  bloq 
   
   
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
         [w, ww] =  lan_model_stat(FT(:,:,T1:T2,:),ones(size(NGO(~noLAN))),NGO(~noLAN), nGO(~noLAN),nNGO(~noLAN),GOOD(~noLAN),pGOOD(~noLAN),LEFT(~noLAN), cfgM);
         
         LAN.freq.model.t =ww.t;
         LAN.freq.model.p= w;
         LAN.freq.model.b= ww.b;
         LAN.freq.model.r= {'Int' ,'NGO','nGO','nNGO','GOOD','pGOOD','LEFT'};
         LAN.freq.model.freq=LAN.freq.freq;
         LAN.freq.model.time=LAN.freq.time(T1:T2);
         
         LAN.data=[];
         
         save ([ PATH_MAT 'LAN_' TAREA '_MODEL_NGO_nGO_nNGP_G_pG_Lf' ],'LAN', '-v7.3')
         save(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG/' DATE '_pro/LAN_' TAREA '_MODEL_NGO_nGO_nNGP_G_pG_Lf'] , 'LAN', '-v7.3')

         cfgM=[];
         cfgM.type = 'glm';
         cfgM.ops = ['pre(glm:_M1_' SU '_' TAREA ' + Noise )'];
         [w, ww] =  lan_model_stat(FT(:,:,T1:T2,:),ones(size(NGO(~noLAN))),NGO(~noLAN), nGO(~noLAN),nNGO(~noLAN),GOOD(~noLAN),pGOOD(~noLAN),LEFT(~noLAN),NoiseLevel(~noLAN),LowSignal(~noLAN),cfgM);         
         LAN.freq.model.t =ww.t;
         LAN.freq.model.p= w;
         LAN.freq.model.b= ww.b;
         LAN.freq.model.r= {'Int' ,'NGO','nGO','nNGO','GOOD','pGOOD','LEFT','noiseNL', 'noiseLS'};
         LAN.freq.model.freq=LAN.freq.freq;
         LAN.freq.model.time=LAN.freq.time(T1:T2);
         
         LAN.data=[];
    
         save ([ PATH_MAT 'LAN_' TAREA '_MODEL_NGO_nGO_nNGP_G_pG_Lf_noise2' ],'LAN', '-v7.3')
         save(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG/' DATE '_pro/LAN_' TAREA '_MODEL_NGO_nGO_nNGP_G_pG_Lf_noise2' ] , 'LAN', '-v7.3')
        
         

end
end


% %% for checking data de data is OK
if 0 
figure
pcolor2(LAN.freq.model.time, LAN.freq.model.freq, squeeze(LAN.freq.model.t{6}(:,59,:))), shading flat
end


