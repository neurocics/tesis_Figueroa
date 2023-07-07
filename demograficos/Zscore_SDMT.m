% Datos de SDMT

puntajeObservado = 57;  % Puntaje observado
media = 57.41  % Media
desviacionEstandar = 12.66;  % Desviación estándar

% Cálculo del puntaje Z
puntajeZ = (puntajeObservado - media) / desviacionEstandar;

% Mostrar el resultado
fprintf('El puntaje Z es: %.2f\n', puntajeZ);