 % Fechas de diagnóstico MS
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

%%


SU={'B_FE_24011989' % 'RL GN WM' 
} ;

GR={'controles'};


D = ls ( ['/Volumes/Alehermosa/TESIS/' GR{1} '/' SU{1} '/EEG/*GN.eeg'])

fechaEEG = datetime([ D(end-11:end-8)  '-' D(end-13:end-12) '-'  D(end-15:end-14)]);

fechaNAC = datetime([ D(end-20:end-17)  '-' D(end-22:end-21) '-'  D(end-24:end-23)]);

% calcular edad de participantes en el momento del EEG MS

edad = years(fechaEEG - datetime(fechaNAC, 'InputFormat', 'yyyy-MM-dd'));

fprintf('Tu edad en el día del experimento es: %.2f años\n', edad);

