% read behave

clear

% lista suejto 

controles={
'B_FE_24011989'
'CAHG_27061988'
'CEAS_23071992'
'GRCN_21111979'
'KXMS_13061971'
'MAAO_15081988'
'MEFR_14061991'
'SNCS_22121989'
'YDCL_11021994'
}';

pacientes={
'S_WL_25051979'
'CAPG_27061969'
'GAMY_03071978'
'KSDG_20051983'
'MJDD_01101984'
'MMPM_06011977'
'CAMS_24101973'
'CFOI_23101987'
'FLFS_08081990'
'K_BA_22111991'
'MPRM_13061962'
'PMRF_18051979'
'TPAO_16051980'
}';

SU_all={controles,pacientes};
GR_all={'controles','pacientes'};
%%
warning off
for g = 1:2
    
    Sub=SU_all{g};
    GR=GR_all{g};
for s = 1:numel(Sub)
    % s =1
    
    SU = Sub{s};
    PATH_D = '/Volumes/Alehermosa/TESIS/%G/%S/LOG/';
    PATH_D = strrep( strrep(PATH_D,'%G',GR) , '%S' , SU);
    
    cd([ PATH_D])
    disp(Sub{s})
    
    file = ls_lan ('*gonogo_lum_sin_cali.log'); 
 
    file = file{~isemptycell(file)};
    
    cfg=[];
    cfg.filename = file;
    cfg.type        = 'presentation';
        cfg.est         = [ 10 21];
        cfg.resp        = [1 2 3 4] ;
        cfg.rw          = [10000] ;      %   (ms)  
        
        RT = rt_read(cfg);
        RT.OTHER.SU=repmat(Sub(s),size(RT.est));
        RT.OTHER.GR=repmat({GR},size(RT.est));
        
        if isempty(RT.OTHER.SU), error(); end
        RT.OTHER = rmfield(rmfield(rmfield(rmfield(RT.OTHER    ,'misport_time_error'),'misport_time'),'misport_code'),'misport_port');
        
        if s==1 && g==1
            TRT=RT;
        else
            TRT=rt_merge(TRT,RT);
        end
    
        % RL
        
   
    file = ls_lan ('*RvL_eye_sin_cali.log'); 

    file = file{~isemptycell(file)};
        
        
     cfg             = [];
        cfg.filename    = file;% '../LOG/OAPC_01041962_23062021_WM-Sternberg_EEG__no_eye_tracking_lum.log';% 
        cfg.type        = 'presentation';
        cfg.est         = [30 31];
        cfg.resp        = [1 2 3 4] ;
        cfg.rw          = [10000] ;      %   (ms)    
        RT_lr = rt_read(cfg); 
        
        
        RT_lr.OTHER.SU=repmat(Sub(s),size(RT_lr.est));
        RT_lr.OTHER.GR=repmat({GR},size(RT_lr.est));
        
        
        if isempty(RT_lr.OTHER.SU), error(); end
        
        file = ls_lan ('*RL_output_RvL_file.txt');
        file = file{~isemptycell(file)};
        OUT_FILE =readtable(file,...
         'Delimiter','\t','ReadVariableNames',true);  

        RT_lr.OTHER.out_latency = OUT_FILE.ttp;    
        RT_lr.OTHER.out_rt = OUT_FILE.ttr - OUT_FILE.ttp;  
        RT_lr.OTHER.out_Bd1 = OUT_FILE.Bd1;	
        RT_lr.OTHER.out_Bd2 = OUT_FILE.Bd2;	        
        RT_lr.OTHER.out_buena = OUT_FILE.buena;	       
        RT_lr.OTHER.out_resp = OUT_FILE.Resp;
        RT_lr.OTHER.out_F = OUT_FILE.F;
        RT_lr.OTHER.out_prob = OUT_FILE.prob;
        RT_lr.OTHER.out_nt = OUT_FILE.nt;
        
        
        RT_lr.OTHER = rmfield(rmfield(rmfield(rmfield(RT_lr.OTHER    ,'misport_time_error'),'misport_time'),'misport_code'),'misport_port');
        
        if s==1 && g==1
            TRT_lr=RT_lr;
        else
            TRT_lr=rt_merge(TRT_lr,RT_lr);
        end       
       
       % WM
       
    file = ls_lan ('*Sternberg_EEG__eye_tracking_lum_sin_cali.log'); 

    file = file{~isemptycell(file)};
        
        
        cfg             = [];
        cfg.filename    = file;% '../LOG/OAPC_01041962_23062021_WM-Sternberg_EEG__no_eye_tracking_lum.log';% 
        cfg.type        = 'presentation';
        cfg.est         = [20 21 40 41 60 61];
        cfg.resp        = [1 2] ;
        cfg.rw          = [10000] ;      %   (ms)    
        RT_wm= rt_read(cfg);
        RT_wm.OTHER.SU=repmat(Sub(s),size(RT_wm.est));
        RT_wm.OTHER.GR=repmat({GR},size(RT_wm.est));
        
           % Correct responses 
   % Succefull Memory performance 
   nB1=90;
   nB2=90;
   if length(RT_wm.est)<180 && length(RT_wm.est)>90
       %error('fixme')
       %continue
       nB1=length(RT_wm.est)-90;
   elseif length(RT_wm.est)==90
       nB1=90;
       nB2=0;
   elseif length(RT_wm.est)<90
       %error('fixme')
       %continue 
       nB1=length(RT_wm.est);
       nB2=0;
   end
   % check responces per  bloq 
   B1 =  [ ones([1 nB1])  zeros([ 1 nB2])  ];
   b1r20 = RT_wm.resp(RT_wm.est==20 & B1 );
   b1r21 = RT_wm.resp(RT_wm.est==21 & B1 );
   b1r40 = RT_wm.resp(RT_wm.est==40 & B1 );
   b1r41 = RT_wm.resp(RT_wm.est==41 & B1 );
   b1r60 = RT_wm.resp(RT_wm.est==60 & B1 );
   b1r61 = RT_wm.resp(RT_wm.est==61 & B1 );
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
   b2r20 = RT_wm.resp(RT_wm.est==20 & B2 );
   b2r21 = RT_wm.resp(RT_wm.est==21 & B2);
   b2r40 = RT_wm.resp(RT_wm.est==40 & B2 );
   b2r41 = RT_wm.resp(RT_wm.est==41 & B2);
   b2r60 = RT_wm.resp(RT_wm.est==60 & B2 );
   b2r61 = RT_wm.resp(RT_wm.est==61 & B2);
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
        
    SMP_b1_1 = (RT_wm.est==21|RT_wm.est==41|RT_wm.est==61) & (B1c==RT_wm.resp); 
    SMP_b1_0 = (RT_wm.est==20|RT_wm.est==40|RT_wm.est==60) & (B1c0==RT_wm.resp); 
    SMP_b1 = SMP_b1_1 +SMP_b1_0;
    
    SMP_b2_1 = (RT_wm.est==21|RT_wm.est==41|RT_wm.est==61) & (B2c==RT_wm.resp); 
    SMP_b2_0 = (RT_wm.est==20|RT_wm.est==40|RT_wm.est==60) & (B2c0==RT_wm.resp); 
    SMP_b2 = SMP_b2_1 + SMP_b2_0;
    
    RT_wm.OTHER.SMP  = [SMP_b1(1:nB1)  SMP_b2(nB1+1:end)];  
   fprintf(['------------------------\n']);
    fprintf(['TOTAL\t\t '   num2str(mean( RT_wm.OTHER.SMP))            '\t\t\n']);
  fprintf(['------------------------\n']);
   
   
        
        if isempty(RT_wm.OTHER.SU), error(); end
        try
        RT_wm.OTHER = rmfield(rmfield(rmfield(rmfield(RT_wm.OTHER    ,'misport_time_error'),'misport_time'),'misport_code'),'misport_port');
        end
        if s==1 && g==1
            TRT_wm=RT_wm;
        else
 
            TRT_wm=rt_merge(TRT_wm,RT_wm);
        end         

        
        
       
       
end

end
%%

cd  /Volumes/Alehermosa/TESIS
COR.RT = TRT;
cfg.filename='COR_gn.txt';
COR2tableR(COR,cfg)

COR.RT = TRT_wm;
cfg.filename='COR_wm.txt';
COR2tableR(COR,cfg)

COR.RT = TRT_lr;
cfg.filename='COR_rl.txt';
COR2tableR(COR,cfg)
