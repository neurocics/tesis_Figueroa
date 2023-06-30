% 
% Model RL from 
% Time-frequiency amplitude 
%


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


for nG = 1:2
    GR= GRS{nG};
    Sujetos_g = Sujetos{nG};
for  nS =Sujetos_g %%
    
for nT = {'RL'}
   disp(['Sujeto: ' nS{1} '(' GR ')'  ' Tarea: '  nT{1} ])
   
   
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
    elseif ~OVERWRITE &&  exist ([ PATH_MAT 'LAN_' TAREA '_MODEL_DM_r_Uc_LD_Ri.mat' ], 'file')
      fprintf(['freq OK ' SU ' \n'])  
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
  
    % periodos diferentes para las regresiones 
%     switch SU
%         case 'IABM_03061982'     
%         otherwise 
%         Feedback = (LAN.RT.est==6)|(LAN.RT.est==7);
%         DM = (LAN.RT.est==4)|(LAN.RT.est==5);
%     end
   Feedback = (LAN.RT.est==6)|(LAN.RT.est==7)|(LAN.RT.est==30)|(LAN.RT.est==31);
   DM = (LAN.RT.est==4)|(LAN.RT.est==5)|(LAN.RT.est==10);
   if any((DM+Feedback)~=1)
        error('Data no calsa!')
    end  
    % sacar trails eliminados 
    DMft = DM(~noFT);
    DM(noLAN)= 0;
    Feedbackft = Feedback(~noFT);
    Feedback(noLAN)= 0;
    
    
    % time
    T1 = find_approx(LAN.freq.time, -0.5);
    T2 = find_approx(LAN.freq.time, 2);
   
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
    
    
    %%% modelo conductual por sujeto 
    FILE_txt =ls_lan([ PATH_LOG '*_RL_output_RvL_file.txt']);
    Tabla =readtable(FILE_txt{1},...
    'Delimiter','\t','ReadVariableNames',true);


    %Searching misis firts trials  
    if length(LAN.RT.OTHER.out_Bd1) < 2*length(Tabla.Bd1)
        dif = 2*length(Tabla.Bd1) - length(LAN.RT.OTHER.out_Bd1);
        Bd1 = repelem(Tabla.Bd1,2);
        if any(LAN.RT.OTHER.out_Bd1 ~=Bd1(dif+1:end)')
            error(' no encontre la diferencia de ensayos ')
        else
            dif2=dif;
            dif= dif/2;
        end
    elseif  length(LAN.RT.OTHER.out_Bd1) == 2*length(Tabla.Bd1)
            dif2=0;
            dif= 0;
    else
         error(' no encontre la diferencia de ensayos ')     
    end
    
    
                if 1 %ifjags
                nchains  = 3; % How Many Chains?
                nburnin  = 1000; % How Many Burn-in Samples?
                nsamples = 1000;  % How Many Recorded Samples?


                firts = find(([1 diff(Tabla.prob)'] ~= 0 )|([1 diff(Tabla.buena)'] ~= 0 ));
                firstT=zeros(size(Tabla.nt));
                firstT(firts(1:2:end))=1;

                % Create a single structure that has the data for all observed JAGS nodes
                DJ.nT = numel(Tabla.nt);
                DJ.feedback = Tabla.F;
                DJ.y = single(Tabla.Resp==1);
                DJ.Bd1 = Tabla.Bd1;
                DJ.Bd2 = Tabla.Bd2;
                DJ.firstT = firstT;

                % Set initial values for latent variable in each chain
                init0=[];
                init0(1).beta0=0;init0(2).beta0=-1;init0(3).beta0=1;
                init0(1).beta1=0;init0(2).beta1=1;init0(3).beta1=1; 
                init0(1).alpha_r=0.5;init0(2).alpha_r=0.7;init0(3).alpha_r=0.5; 
                init0(1).alpha_a=0.2;init0(2).alpha_a=0.4;init0(3).alpha_a=0.3;
                init0(1).rew=0.2;init0(2).rew=1;init0(3).rew=0.5;
                init0(1).gamma=0.2;init0(2).gamma=1;init0(3).gamma=0.5;

                % fitting
                %  U[i] <- V.1[i]  * Pi.1[i]  - V.2[i] * Pi.2[i]
                tic
                doparallel = 0;
                fprintf( 'Running JAGS\n' );
                [samples, stats ] = matjags( ...
                        DJ, ...
                        fullfile('/Users/alejandra/Desktop/DOCTORADO/2021/TESIS /git_neurocovid/tesis_Figueroa/EEG-EYE', 'model_rL_learn_1a.txt'), ...
                        init0, ...
                        'doparallel' , doparallel, ...
                        'nchains', nchains,...
                        'nburnin', nburnin,...
                        'nsamples', nsamples, ...
                        'thin', 10, ...
                        'monitorparams', {'beta0','beta1','alpha_r','alpha_a','rew','delta_a','delta_r','U','V.1','Pi.1','V.2','Pi.2'  }, ...
                        'savejagsoutput' , 0 , ...
                        'verbosity' , 1 , ...
                        'cleanup' , 0  );
                toc

                RT.OTHER.Resp = repelem(Tabla.Resp,2)';
                
                RT.OTHER.LearningSignal = (stats.mean.delta_a(2:end).*Tabla.F')+(stats.mean.delta_r(2:end).*(1-Tabla.F'));
                RT.OTHER.LearningSignal_dm = [0 RT.OTHER.LearningSignal(1:end-1)];
                RT.OTHER.LearningSignal_dm_c = [0 RT.OTHER.LearningSignal(1:end-1)].* (-1*(Tabla.Resp==2)' + 1*(Tabla.Resp==1)');

                RT.OTHER.LearningSignal = repelem(RT.OTHER.LearningSignal,2);
                RT.OTHER.LearningSignal_dm = abs(repelem(RT.OTHER.LearningSignal_dm,2));
                RT.OTHER.LearningSignal_dm_c = repelem(RT.OTHER.LearningSignal_dm_c,2);

                RT.OTHER.PredictionError= (stats.mean.delta_a(2:end).*Tabla.F'/stats.mean.alpha_a)+(stats.mean.delta_r(2:end).*(1-Tabla.F')/stats.mean.alpha_a);             
                RT.OTHER.PredictionError = RT.OTHER.PredictionError .* (-1*(Tabla.Resp==2) + 1*(Tabla.Resp==1))';
                RT.OTHER.PredictionError = repelem(RT.OTHER.PredictionError,2);

                RT.OTHER.Utility_c = stats.mean.U .* (-1*(Tabla.Resp==2)' + 1*(Tabla.Resp==1)');
                RT.OTHER.Utility_c = repelem(RT.OTHER.Utility_c,2);
                
                RT.OTHER.Utility_z = normal_z(RT.OTHER.Utility_c);
                RT.OTHER.LearningSignal_dm_z = normal_z(RT.OTHER.LearningSignal_dm);
               
                % truncate extreme values
                RT.OTHER.Utility_z(RT.OTHER.Utility_z>3.5)  = 3.6;
                RT.OTHER.Utility_z(RT.OTHER.Utility_z<-3.5)  = -3.6;
                RT.OTHER.LearningSignal_dm_z(RT.OTHER.LearningSignal_dm_z>3.5)  = 3.6;
                RT.OTHER.LearningSignal_dm_z(RT.OTHER.LearningSignal_dm_z<-3.5)  = -3.6;
                
                
                RT.OTHER.V1 = repelem( stats.mean.V_1 , 2);
                RT.OTHER.V2 = repelem( stats.mean.V_2 , 2);
                RT.OTHER.Vc = (RT.OTHER.V1-RT.OTHER.V2) .* (-1*(RT.OTHER.Resp==2) + 1*(RT.OTHER.Resp==1));
                
                
                RT.OTHER.Pi1 = repelem( stats.mean.Pi_1 , 2);
                RT.OTHER.Pi2 = repelem( stats.mean.Pi_2 , 2);
                RT.OTHER.Pic = (RT.OTHER.Pi1-RT.OTHER.Pi2) .* (-1*(RT.OTHER.Resp==2) + 1*(RT.OTHER.Resp==1));
                for R=fields(RT.OTHER)' 
                    eval(['RT.OTHER.' R{1} ' = RT.OTHER.' R{1} '(' num2str(dif2+1) ':end) ;' ])
                end
                
                end
    
                save ([ PATH_LOG 'MODELO_' TAREA '_1a_EEG_' ], 'RT', 'stats', 'samples');
    
      Resp_IZQ = (LAN.RT.OTHER.out_resp==1);
    
    
    
    
    
    
    %%% fitting model 
    
         cfgM=[];
         cfgM.type = 'glm';
         cfgM.ops = ['pre(glm:_M1_' SU  ' _ ' TAREA '  )'];
         [w, ww] =  lan_model_stat(FT(:,:,T1:T2,DMft),RT.OTHER.Utility_c(DM), RT.OTHER.LearningSignal_dm(DM),Resp_IZQ(DM) ,cfgM);
         
         LAN.freq.model.t =ww.t;
         LAN.freq.model.p= w;
         LAN.freq.model.b= ww.b;
         LAN.freq.model.r= {'Int' ,'Utility_c','LS','R_izq'};
         LAN.freq.model.freq=LAN.freq.freq;
         LAN.freq.model.time=LAN.freq.time(T1:T2);
         
         LAN.data=[];
  
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
         
         
         
         
         
         save ([ PATH_MAT 'LAN_' TAREA '_MODEL_DM_r_Uc_LD_Ri' ],'LAN', '-v7.3')
         % some problems with drive 
         try
                  %save(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG/' DATE '_pro/LAN_' TAREA '_MODEL_DM_r_Uc_LD_Ri' ] , 'LAN', '-v7.3')
         catch
                 pause(30)
                 try
                 %save(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG/' DATE '_pro/LAN_' TAREA '_MODEL_DM_r_Uc_LD_Ri' ] , 'LAN', '-v7.3')
                 end
         end
         
         cfgM=[];
         cfgM.type = 'glm';
         cfgM.ops = ['pre(glm:_M1_' SU ' _ ' TAREA ' + Noise )'];
         [w, ww] =  lan_model_stat(FT(:,:,T1:T2,DMft),RT.OTHER.Utility_c(DM), RT.OTHER.LearningSignal_dm(DM),Resp_IZQ(DM),NoiseLevel(DM),LowSignal(DM),cfgM);
         
         LAN.freq.model.t =ww.t;
         LAN.freq.model.p= w;
         LAN.freq.model.b= ww.b;
         LAN.freq.model.r= {'Int' ,'Utility_c','LS','R_izq','noiseNL', 'noiseLS'};
         LAN.freq.model.freq=LAN.freq.freq;
         LAN.freq.model.time=LAN.freq.time(T1:T2);
         
         LAN.data=[];
    
         save ([ PATH_MAT 'LAN_' TAREA '_MODEL_DM_r_Uc_LD_Ri_noise2' ],'LAN', '-v7.3')
         
         % some problems with drive 
         try
         %save(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG/' DATE '_pro/LAN_' TAREA '_MODEL_DM_r_Uc_LD_Ri_noise2' ] , 'LAN', '-v7.3')
         catch
             pause(30)
             try
                 %save(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG/' DATE '_pro/LAN_' TAREA '_MODEL_DM_r_Uc_LD_Ri_noise2' ] , 'LAN', '-v7.3')
             end
         end
         

    






end
end
end

% %% for checking data de data is OK
if 0 
figure
pcolor2(LAN.freq.model.time, LAN.freq.model.freq, squeeze(LAN.freq.model.t{6}(:,59,:))), shading flat
end

