% primera lectura de analisis de datos de modelamiento y conducta

clear

SUBS = { 'MJDD_01101984'}% sin datos 
    %'UENM_22121979'}
 





for s = 1:length(SUBS)

SU=SUBS{s};% SU='JFFA_04101976'
GR='';

PATH_D = '/Volumes/Alehermosa/TESIS/%G/%S/LOG/';
PATH_D = strrep( strrep(PATH_D,'%G',GR) , '%S' , SU);

%PATH_Google = '/Volumes/GoogleDrive/Mi unidad/DATA/PARTICIPANTES/%S/LOG' ;
%PATH_Google = strrep( strrep(PATH_Google,'%G',GR) , '%S' , SU);




if 1 % exist(PATH_D,'dir')==7


%cd(PATH_D )
cd(PATH_Google)
data_log = ls_lan(['*gonogo_lum*.log' ]);
data_log(logical(fun_in_cell(data_log,'isempty(@)')))=[];
data_log(logical(fun_in_cell(data_log,'contains(@,''training'')')))=[];
data_log = data_log{1};%(1:end-1);

        % .log
        cfg             = [];
        cfg.filename    = data_log;% '../LOG/OAPC_01041962_23062021_WM-Sternberg_EEG__no_eye_tracking_lum.log';% 
        cfg.type        = 'presentation';
        cfg.est         = [ 10 21];
        cfg.resp        = [1 2] ;
        cfg.rw          = [10000] ;      %   (ms) 
        warning off
        
RT_log = rt_read(cfg);
warning on 

RT_log.OTHER.Sub = repmat({SU},size(RT_log.est));
seqG = zeros(size(RT_log.est));

for t = 1:length(RT_log.est)
    if t==1 && RT_log.est(t)==10
                seqG(t) = 1 ;
    elseif t>1 && RT_log.est(t)==10 
                seqG(t) = 1 +  seqG(t-1);
    elseif  RT_log.est(t)==21
                seqG(t)= 0;
    end
end

RT_log.OTHER.seqG = seqG;

if s == 1
    RT_all = RT_log;
else
    RT_all = rt_merge(RT_all, RT_log);
end


% RL task 
cd(PATH_Google)
data_txt = ls(['*RL_output_RvL_file.txt' ]);
data_txt = data_txt(1:end-1);

tabla = readtable(data_txt);

RT_txt  = [];
RT_txt.est = tabla.nt';
RT_txt.rt = tabla.ttr' - tabla.ttp';
RT_txt.laten = tabla.ttp';
RT_txt.resp = tabla.Resp';
RT_txt = rt_check(RT_txt);

RT_txt.OTHER.Sub = repmat({SU},size(RT_txt.est));
RT_txt.OTHER.nt = tabla.nt';
RT_txt.OTHER.buena = tabla.buena';
RT_txt.OTHER.Bd1 = tabla.Bd1';
RT_txt.OTHER.Bd2 = tabla.Bd2';
RT_txt.OTHER.Resp = tabla.Resp';
RT_txt.OTHER.F = tabla.F';
RT_txt.OTHER.prob = tabla.prob';
RT_txt.OTHER.ttp = tabla.ttp';
RT_txt.OTHER.ttr = tabla.ttr';
RT_txt.OTHER.ttF = tabla.ttF';

if s == 1
    RT_rl = RT_txt;
else
    RT_rl = rt_merge(RT_rl, RT_txt);
end



else
   disp( ['no hay datos para: ' SU '<----------']) 
    
end

clear RT_log

end

%%
fprintf('\ndatos cargados\n')

cd('/Volumes/GoogleDrive/Mi unidad/DOCTORADO/2021/UI/DATA/')

COR.RT = RT_all;
%

% write data into .txt for R 
cfg=[];
cfg.filename = 'data_GN_neuroCOVID.txt' ; %       
cfg.format  =  'txt';% 'mat'
cfg.delimiter = '\t';
fprintf('\ngrabando Gonogo ...\n')
COR2tableR(COR,cfg)

COR= [];
COR.RT = RT_rl;
%

% write data into .txt for R 
cfg=[];
cfg.filename = 'data_RL_eeg_neuroCOVID.txt' ; %       
cfg.format  =  'txt';% 'mat'
cfg.delimiter = '\t';
fprintf('\ngrabando RL ... \n')
COR2tableR(COR,cfg)

fprintf('\nListo ama .... :) \n')