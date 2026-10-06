%% A4 - HMM Gaussiano: cada entrada frente a todas las salidas
clear; clc; close all;
cfg.root = 'F:\Datos de investigacion\Tomate\Datos sincronos\Version 3';
cfg.regime = menu('Seleccione la base','1 - Dia','2 - Noche','3 - Madrugada');
cfg.kRange = 2:8; cfg.replicates = 2; cfg.maxIter = 100; cfg.tol = 1e-5; cfg.randomSeed = 42;
cfg.nPCsMax = 6; cfg.pcaVariance = 90; cfg.minVariance = 1e-4; cfg.transitionPseudoCount = 1e-3; cfg.saveFigures = true;
A4 = A4_Funcion_HMM(cfg);
