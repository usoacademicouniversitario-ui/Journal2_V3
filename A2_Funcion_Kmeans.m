function Result = A2_Funcion_Kmeans(cfg)
% A2: K-means para cada ENTRADA + TODAS las SALIDAS. No existe clustering de salidas solas.
rng(cfg.randomSeed,'twister');D=local_load(cfg);out=fullfile(cfg.root,'Resultados_A',D.regime,'A2_Kmeans');if ~isfolder(out),mkdir(out);end
summary=table();Result=struct();opts=statset('MaxIter',cfg.maxIter,'Display','off');
for i=1:numel(D.inputs)
 names=[D.inputs(i),D.outputs];Xi=D.X(:,[i,numel(D.inputs)+(1:numel(D.outputs))]);mask=all(isfinite(Xi),2);Xc=Xi(mask,:);if size(Xc,1)<100,continue;end
 [Z,~,~,~,keep]=local_scale(Xc,names);namesK=names(keep);if ~keep(1),warning('Entrada %s constante: A2 omitido.',D.inputs{i});continue;end
 sil=NaN(numel(cfg.kRange),1);wcss=sil;models=cell(numel(cfg.kRange),1);
 for q=1:numel(cfg.kRange),k=cfg.kRange(q);[idx,C,sumd]=kmeans(Z,k,'Distance','sqeuclidean','Replicates',cfg.replicates,'Options',opts,'Start','plus');sv=silhouette(Z,idx,'sqeuclidean');sil(q)=mean(sv);wcss(q)=sum(sumd);models{q}=struct('idx',idx,'C',C,'sv',sv);end
 [~,b]=max(sil);bestK=cfg.kRange(b);M=models{b};tag=char(local_safe(D.inputs{i}));folder=fullfile(out,tag);if ~isfolder(folder),mkdir(folder);end
 writetable(table(cfg.kRange(:),sil,wcss,'VariableNames',{'K','MeanSilhouette','WCSS'}),fullfile(folder,'K_Selection.csv'));
 A=table(D.dayID(mask),D.rowID(mask),string(D.timeText(mask)),M.idx,M.sv,'VariableNames',{'DayID','RowInDay','Time','Cluster','Silhouette'});writetable(A,fullfile(folder,'Assignments.csv'));
 R=struct('input',D.inputs{i},'outputs',{D.outputs},'names',{namesK},'bestK',bestK,'idx',M.idx,'centroids',M.C,'silhouette',M.sv,'completeMask',mask);save(fullfile(folder,'A2_Kmeans_Result.mat'),'R','-v7.3');
 summary=[summary;table(string(D.inputs{i}),size(Xc,1),bestK,sil(b),'VariableNames',{'Input','CompleteRows','BestK','MeanSilhouette'})];
 if cfg.saveFigures,f=figure('Visible','off','Color','w');plot(cfg.kRange,sil,'-o');xlabel('K');ylabel('Mean silhouette');grid on;exportgraphics(f,fullfile(folder,'Kmeans_K.png'),'Resolution',180);close(f);end
end
writetable(summary,fullfile(out,'A2_Summary.csv'));Result.summary=summary;Result.regime=D.regime;Result.outputFolder=out;save(fullfile(out,'A2_Resultados.mat'),'Result','-v7.3');fprintf('\nA2 K-means terminado | %s | %d entradas\n%s\n',D.regime,height(summary),out);
end

function D = local_load(cfg)
labels={'Dia','Noche','Madrugada'}; if ~isscalar(cfg.regime)||~ismember(cfg.regime,1:3),error('Seleccion invalida.');end
reg=labels{cfg.regime}; file=fullfile(cfg.root,sprintf('AGC4_Trigger_%s_69_Jornadas.xlsx',reg)); if ~isfile(file),error('No se encontro: %s',file);end
sheets=sheetnames(file); if numel(sheets)~=69,warning('Se esperaban 69 hojas y se encontraron %d.',numel(sheets));end
inputs={'weather/air_temperature.outside','weather/relative_humidity.outside','weather/radiation_global','weather/wind_speed','weather/rain_state','weather/par.outside','weather/heat_emission','compartment/heating_lower_circuit/pipe_temperature','compartment/screen_energy/screen_position','compartment/screen_blackout/screen_position','compartment/water_supply/water_flow_duration','compartment/window_position_lee_side','compartment/window_position_wind_side','compartment/co2_actuation_state','compartment/lamps_activation_percentage'};
outputs={'compartment/water_drain/water_volume','compartment/water_drain/ec','compartment/air_temperature','compartment/relative_humidity','compartment/par','compartment/co2_concentration','compartment/substrate/relative_permittivity','compartment/substrate/bulk_ec','compartment/substrate/substrate_temperature','compartment/substrate/relative_permittivity.1','compartment/substrate/bulk_ec.1','compartment/substrate/substrate_temperature.1','compartment/substrate/relative_permittivity.2','compartment/substrate/bulk_ec.2','compartment/substrate/substrate_temperature.2','compartment/leaf_temperature','compartment/leaf_temperature.1'};
X=[];dayID=[];rowID=[];timeText={};
for s=1:numel(sheets)
 T=readtable(file,'Sheet',sheets{s},'VariableNamingRule','preserve'); need=[{'hora'},inputs,outputs];miss=setdiff(need,T.Properties.VariableNames);if ~isempty(miss),error('Hoja %s: faltan columnas: %s',sheets{s},strjoin(miss,', '));end
 M=NaN(height(T),numel(inputs)+numel(outputs));for j=1:numel(inputs),M(:,j)=double(T.(inputs{j}));end;for j=1:numel(outputs),M(:,numel(inputs)+j)=double(T.(outputs{j}));end
 ir=find(strcmp(inputs,'compartment/water_supply/water_flow_duration'));d=[NaN;diff(M(:,ir))];d(d<0)=NaN;M(:,ir)=d;
 X=[X;M];dayID=[dayID;repmat(s,height(T),1)];rowID=[rowID;(1:height(T))'];h=string(T.hora);timeText=[timeText;cellstr(h)];
end
inputs{strcmp(inputs,'compartment/water_supply/water_flow_duration')}='IrrigationIncrement';
D=struct('file',file,'regime',reg,'X',X,'inputs',{inputs},'outputs',{outputs},'dayID',dayID,'rowID',rowID,'timeText',{timeText});
end
function [Z,med,scale,method,keep]=local_scale(X,names)
keep=true(1,size(X,2));med=NaN(1,size(X,2));scale=NaN(1,size(X,2));method=strings(1,size(X,2));
for j=1:size(X,2),x=X(:,j);med(j)=median(x);q=prctile(x,[25 75]);iq=q(2)-q(1);sd=std(x);if isfinite(iq)&&iq>0,scale(j)=iq;method(j)="IQR";elseif isfinite(sd)&&sd>0,scale(j)=sd;method(j)="STD_FALLBACK";else,keep(j)=false;method(j)="CONSTANT_SKIPPED";warning('Variable constante omitida: %s',names{j});end,end
Z=(X(:,keep)-med(keep))./scale(keep);
end
function s=local_safe(s),s=regexprep(s,'[^A-Za-z0-9]+','_');s=regexprep(s,'^_+|_+$','');if strlength(s)>60,s=extractBefore(s,61);end,end
