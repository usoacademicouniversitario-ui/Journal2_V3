%% A2 - K-means: cada entrada frente a todas las salidas
clear; clc; close all;
cfg.root = 'F:\Datos de investigacion\Tomate\Datos sincronos\Version 3';
cfg.regime = menu('Seleccione la base','1 - Dia','2 - Noche','3 - Madrugada');
cfg.kRange = 2:8; cfg.replicates = 20; cfg.maxIter = 1000; cfg.randomSeed = 42; cfg.saveFigures = true;
A2 = A2_Funcion_Kmeans(cfg);
