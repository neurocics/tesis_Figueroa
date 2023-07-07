 % Calcular tiempo  de diagnóstico MS al día del experimento (EEG)
 
fechasDiagnostico = ["2019-01-01", "2015-01-01", "2011-01-01", "2013-01-01", "2003-01-01", "2013-01-01", "2012-01-01", "2013-01-01", "2016-01-01", "2018-01-01", "2022-01-01", "2016-01-01"]; % Año diagnóstico

% Obtener la fecha actual

fechaActual = datetime('2023-01-01');

% Calcular los años de diagnostico para cada paciente
edades = years(fechaActual - datetime(fechasDiagnostico, 'InputFormat', 'yyyy-MM-dd'));

% Mostrar los años de diagnostico para cada paciente
disp("años diagnostico de MS:");
for i = 1:length(edades)
    disp(['Paciente', num2str(i), ': ', num2str(edades(i)), ' años']);
end

% save en carpeta demografica 
save '/Volumes/Alehermosa/TESIS/demo'

%% Calcular edad de participantes al día del experimento (EEG)


SU= {{
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
'SFVB_26031991'
'N_LS_15111975'
'FSDT_10101990'
'DITR_30081986'
'J_AA_31071995'
'MALV_18031978'
'TPAO_16051980'
'PMRF_18051979'
'MPRM_13061962'

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
'PIHV_23091978'
'JFFA_04101976'
'PGTR_30081986'
'AAAS_13071969'
'MJMM_18101984'
'PSAV_10041994'
}'} ;


GRS={'pacientes','controles'};

OVERWRITE=false;

%%
% Abrir el archivo en modo escritura (si el archivo no existe, se creará)
archivo = fopen('/Volumes/Alehermosa/TESIS/datos_demograficos.txt', 'w');

% Verificar si el archivo se abrió correctamente
if archivo == -1
    error('No se pudo abrir el archivo.');
end

for nG = 1:2
    GR= GRS{nG};
    SU_g = SU{nG};
for nS = 1:numel(SU_g)
    % for nS = Su_g
    
D = ls (['/Volumes/Alehermosa/TESIS/' GR '/' SU_g{nS}  '/EEG/*GN.eeg']);
%D = ls (['/Volumes/Alehermosa/TESIS/' GR '/' nS{1}  '/EEG/*GN.eeg']);

fechaEEG = datetime([ D(end-11:end-8)  '-' D(end-13:end-12) '-'  D(end-15:end-14)]);

fechaNAC = datetime([ D(end-20:end-17)  '-' D(end-22:end-21) '-'  D(end-24:end-23)]);

% calcular edad de participantes en el momento del EEG MS

edad = years(fechaEEG - datetime(fechaNAC, 'InputFormat', 'yyyy-MM-dd'));

% calcular edad de participantes en el momento del EEG MS
fprintf('El sujeto %s qu es %s  tiene una edad de: %.1f años\n', SU_g{nS},GR, edad);


fprintf(archivo, '%s,%s,%.10f\n', SU_g{nS},GR, edad);

end
end 

fclose(archivo);


