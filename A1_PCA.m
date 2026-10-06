%% A1 - PCA: cada entrada frente a todas las salidas
clear; clc; close all;
cfg.root = 'F:\Datos de investigacion\Tomate\Datos sincronos\Version 3';
cfg.regime = menu('Seleccione la base','1 - Dia','2 - Noche','3 - Madrugada');
cfg.varianceThresholds = [70 80 90 95]; cfg.topLoadings = 10; cfg.saveFigures = true;
A1 = A1_Funcion_PCA(cfg);
