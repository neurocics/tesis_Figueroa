


clear

SUs ={
'ASCS_31121976'  %;'08012022';
'CACB_27091965'  %;'17122021';
 'BCUV_02091982' %;'10122021';
'C_DG_15111995'  %; '30062021';
'CAMO_08111988'  %;'29092021';
 'CARM_14011974'  %;  %19032021_pro
 'CDPR_25031986' %;'29072021';
 'CEBT_30091989' %'03112021';
  'CMCO_17051975' %;'06102021';
 'EERV_01041955' %'25082021'
 'ELFN_09021961' %;'09092021';
 'FABN_19111984'%
 'FJAM_10012001'%
 'GADM_08111983'%
 'GAQJ_21061987'%
 'IABF_16091997'%
 'IABM_03061982'%
 'JAPV_03111995'%
 'JAZV_27081991'%
 'JDRP_06081956'%
 'JFFA_04101976'%
 'JFSO_05061967'%
 'JGRF_24091991'%
 'LAAB_29101982'%
 'LDEV_31011988'%
 'LMRG_21041962'%
 'MAEM_17051957'%
 'MFCM_21011984'%
 'MGRN_22101979'%
 'MLAD_09021957'%
 'MLCB_11041981'%
 'MPMM_22111977'%
 'MTVV_27111974'%
 'OAPC_01041962'%
 'RABG_09091981'%
 'RAPA_19121977'%
 'REOO_10031989'%
 'RFCR_12111978'%
 'RSVT_17091955'%
 'SAAA_10051991'%
 'SAJC_22021972'%
 'SAVS_23061979'%
 %'UENM_22121979'
 'VAPH_11101985'%
};


CF=[];
CF.type =  'Morlet';
CF.fwin  = [5] ;% cycles per time window for the 'type' MultiTaper ('MTaper')
CF.step  = [10] ; %numero de puntos entre ventanas e.g. = 10
CF.foi  = [1:0.25:5 5.5:0.5:12 13:35]; %[f1f2 f3 f4 f5 f6 ... ] ; eje de frecuencias  e.g.: 
CF.keeptrials = 'yes';%  - 'file' - ('file4chan' only for 'Morlet'  type ) 

load chanlocs_65geodesic_HCGSNv10(64)
chanlocs([23,55,61:64]) = [];



for su =1:19%numel(SUs):-1:20
    % su=1
SU = SUs{su}  ;%;

DA = '';

LUGAR = '/Volumes/DBNC_03/neuroCOVID/DATA/';


     if isempty(DA)
        [paso paso2] = unix(['ls  -d '  LUGAR SU '/EEG/*_pro' ]);
         DA = paso2(end-12:end-5);
     end

     
     % re-reference to mean and delete the no-EEG channels
            
      EYE = load([ LUGAR SU '/EEG/'  DA '_pro/LAN_GN_EEG_EYE_Baseline.mat']);

      load([ LUGAR SU '/EEG/'  DA '_pro/LAN_GN_EEG_EYE_Baseline_freqEEG.mat'],'LAN')
     
         % Noise Nuance Regressor     
        %if any (LAN.tag.mat(:)>0)
          LowSignal = mean(LAN.tag.mat>0,1) ;
        %end
        % Global high frequiency 
        H1 =find_approx(LAN.freq.fourierp.freq,15);
        H2 =find_approx(LAN.freq.fourierp.freq,45);
        NoiseLevel = normal_z( squeeze(nanmean(nanmean(LAN.freq.fourierp.data(:,H1:H2,:)        ,1),2))');

        PU = cat(4,EYE.LAN.data{:});
        paso = PU;
        paso(:) = isnan(paso(:)); 
        LOW_EYES_Signal  = squeeze(nanmean(nanmean(paso([69 67],:,:,:)        ,1),2))';
        
        time = timelan(EYE.LAN);
        t1 =find_approx(time,0.95);
        t2 =find_approx(time,1.05);
        PU = normal_z(PU([70 68],t1:t2,:,:));
        PUPIL = squeeze(nanmean(nanmean( PU      ,1),2))';
        
        %clear EYE
        
        FT = LAN.freq.powspctrm;
        noFT= isemptycell(FT);
        noLAN = LAN.accept==0;
        if any(noFT~=noLAN)
            error('Data no calsa!')
        else
            FT =cat(4,FT{:});
        end
        
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
         cfgM.ops = ['pre(glm:_M1_' SU '_ EEG -EYE )'];
         [w, ww] =  lan_model_stat(FT,PUPIL(~noLAN), ones(size(LowSignal(~noLAN))), LowSignal(~noLAN),  NoiseLevel(~noLAN) , LOW_EYES_Signal(~noLAN), cfgM);
         
         LAN.freq.model.t =ww.t;
         LAN.freq.model.p= w;
         LAN.freq.model.b= ww.b;
         LAN.freq.model.r= {'Int' ,'PUPIL','nosieX3'};
         LAN.freq.model.freq=LAN.freq.freq;
         LAN.freq.model.time=LAN.freq.time;
         
         LAN.data=[];
         
         save ([ LUGAR SU '/EEG/'  DA '_pro/LAN_BL_PUPIL_MODEL' ],'LAN', '-v7.3')
         save(['/Volumes/GoogleDrive-103235447575506129142/Mi unidad/DATA/PARTICIPANTES/' SU '/EEG/' DA '_pro/LAN_BL_PUPIL_MODEL_new'] , 'LAN', '-v7.3')

            
            
          
            % save script in the LAN structure
            if isempty(mfilename('fullpath'))
            name = matlab.desktop.editor.getActiveFilename;
            else
            name=[mfilename('fullpath') '.m' ];
            end  
     
 
     
     
end












