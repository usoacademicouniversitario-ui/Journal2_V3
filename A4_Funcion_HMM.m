function Result=A4_Funcion_HMM(cfg)
% A4: HMM Gaussiano diagonal para cada ENTRADA + TODAS las SALIDAS.
% Cada hoja es una secuencia independiente: nunca hay transiciones entre jornadas.
rng(cfg.randomSeed,'twister');D=local_load(cfg);out=fullfile(cfg.root,'Resultados_A',D.regime,'A4_HMM');if ~isfolder(out),mkdir(out);end
summary=table();Result=struct();
for i=1:numel(D.inputs)
 names=[D.inputs(i),D.outputs];Xi=D.X(:,[i,numel(D.inputs)+(1:numel(D.outputs))]);mask=all(isfinite(Xi),2);Xc=Xi(mask,:);days=D.dayID(mask);if size(Xc,1)<100,continue;end
 [Z,~,~,~,keep]=local_scale(Xc,names);namesK=names(keep);if ~keep(1),warning('Entrada %s constante: A4 omitido.',D.inputs{i});continue;end
 Zc=Z-mean(Z,1);[~,S,V]=svd(Zc,'econ');lat=diag(S).^2/(size(Zc,1)-1);ex=100*lat/sum(lat);npc=find(cumsum(ex)>=cfg.pcaVariance,1);if isempty(npc),npc=size(V,2);end;npc=min(cfg.nPCsMax,npc);Y=(Zc*V(:,1:npc))./sqrt(lat(1:npc))';
 [seqS,seqE]=local_sequences(days);BIC=NaN(numel(cfg.kRange),1);LL=BIC;mods=cell(size(BIC));
 for q=1:numel(cfg.kRange),K=cfg.kRange(q);bestLL=-Inf;best=[];for r=1:cfg.replicates,rng(cfg.randomSeed+1000*q+r);try,m=local_fit_hmm(Y,seqS,seqE,K,cfg);if m.logLikelihood>bestLL,bestLL=m.logLikelihood;best=m;end,catch ME,warning('HMM %s K=%d rep=%d: %s',D.inputs{i},K,r,ME.message);end,end;if ~isempty(best),p=(K-1)+K*(K-1)+2*K*npc;LL(q)=bestLL;BIC(q)=-2*bestLL+p*log(size(Y,1));mods{q}=best;end,end
 valid=find(isfinite(BIC));if isempty(valid),continue;end;[~,rr]=min(BIC(valid));b=valid(rr);best=mods{b};K=cfg.kRange(b);[gamma,~]=local_posterior(Y,seqS,seqE,best);state=local_viterbi(Y,seqS,seqE,best);maxP=max(gamma,[],2);
 tag=char(local_safe(D.inputs{i}));folder=fullfile(out,tag);if ~isfolder(folder),mkdir(folder);end
 writetable(table(cfg.kRange(:),LL,BIC,'VariableNames',{'K','LogLikelihood','BIC'}),fullfile(folder,'HMM_Selection.csv'));
 A=table(days,D.rowID(mask),string(D.timeText(mask)),state,maxP,'VariableNames',{'DayID','RowInDay','Time','State','MaxPosterior'});writetable(A,fullfile(folder,'Assignments.csv'));
 R=struct('input',D.inputs{i},'outputs',{D.outputs},'names',{namesK},'bestK',K,'model',best,'state',state,'posterior',gamma,'pcaCoeff',V(:,1:npc),'pcaExplained',ex,'completeMask',mask,'sequenceStart',seqS,'sequenceEnd',seqE);save(fullfile(folder,'A4_HMM_Result.mat'),'R','-v7.3');
 summary=[summary;table(string(D.inputs{i}),size(Xc,1),npc,K,BIC(b),mean(maxP),'VariableNames',{'Input','CompleteRows','PCs','BestK','BIC','MeanMaxPosterior'})];
 if cfg.saveFigures,f=figure('Visible','off','Color','w');plot(cfg.kRange,BIC,'-o');xlabel('K');ylabel('BIC');grid on;exportgraphics(f,fullfile(folder,'HMM_BIC.png'),'Resolution',180);close(f);end
end
writetable(summary,fullfile(out,'A4_Summary.csv'));Result.summary=summary;Result.regime=D.regime;Result.outputFolder=out;save(fullfile(out,'A4_Resultados.mat'),'Result','-v7.3');fprintf('\nA4 HMM terminado | %s | %d entradas\n%s\n',D.regime,height(summary),out);
end

function [ss,se]=local_sequences(day)
chg=[true;diff(day)~=0];ss=find(chg);se=[ss(2:end)-1;numel(day)];
end

function m=local_fit_hmm(Y,ss,se,K,cfg)
[N,D]=size(Y);idx=kmeans(Y,K,'Replicates',3,'Start','plus','MaxIter',300);mu=zeros(K,D);va=zeros(K,D);
for k=1:K,X=Y(idx==k,:);if isempty(X),X=Y(randi(N,max(2,round(N/K)),1),:);end;mu(k,:)=mean(X,1);va(k,:)=max(var(X,1,1),cfg.minVariance);end
pi0=ones(K,1)/K;A=ones(K)/K;prev=-Inf;
for it=1:cfg.maxIter
 logB=local_logB(Y,mu,va);piAcc=zeros(K,1);Aacc=zeros(K);gSum=zeros(K,1);ySum=zeros(K,D);y2Sum=zeros(K,D);ll=0;
 for s=1:numel(ss),a=ss(s);b=se(s);[g,xi,L]=local_fb(logB(a:b,:),pi0,A);ll=ll+L;piAcc=piAcc+g(1,:)';Aacc=Aacc+xi;gs=sum(g,1)';gSum=gSum+gs;ySum=ySum+g'*Y(a:b,:);y2Sum=y2Sum+g'*(Y(a:b,:).^2);end
 pi0=(piAcc+eps)/sum(piAcc+eps);A=Aacc+cfg.transitionPseudoCount;A=A./sum(A,2);mu=ySum./max(gSum,eps);va=y2Sum./max(gSum,eps)-mu.^2;va=max(va,cfg.minVariance);
 if isfinite(prev)&&abs(ll-prev)<=cfg.tol*(1+abs(prev)),break;end;prev=ll;
end
m=struct('pi',pi0,'A',A,'mu',mu,'var',va,'logLikelihood',ll,'iterations',it);
end

function L=local_logB(Y,mu,va)
K=size(mu,1);L=zeros(size(Y,1),K);for k=1:K,L(:,k)=-0.5*sum(log(2*pi*va(k,:))+((Y-mu(k,:)).^2)./va(k,:),2);end
end

function [gamma,xiSum,ll]=local_fb(logB,pi0,A)
[T,K]=size(logB);rowMax=max(logB,[],2);b=exp(logB-rowMax);c=zeros(T,1);alpha=zeros(T,K);
alpha(1,:)=pi0'.*b(1,:);c(1)=sum(alpha(1,:));alpha(1,:)=alpha(1,:)/max(c(1),realmin);
for t=2:T,alpha(t,:)=(alpha(t-1,:)*A).*b(t,:);c(t)=sum(alpha(t,:));alpha(t,:)=alpha(t,:)/max(c(t),realmin);end
beta=ones(T,K);for t=T-1:-1:1,beta(t,:)=(A*(b(t+1,:)'.*beta(t+1,:)'))'/max(c(t+1),realmin);end
gamma=alpha.*beta;gamma=gamma./sum(gamma,2);xiSum=zeros(K);
for t=1:T-1,x=(alpha(t,:)'*(b(t+1,:).*beta(t+1,:))).*A;x=x/max(sum(x,'all'),realmin);xiSum=xiSum+x;end
ll=sum(log(max(c,realmin)))+sum(rowMax);
end

function [G,ll]=local_posterior(Y,ss,se,m)
logB=local_logB(Y,m.mu,m.var);G=zeros(size(Y,1),size(m.A,1));ll=0;for s=1:numel(ss),[g,~,L]=local_fb(logB(ss(s):se(s),:),m.pi,m.A);G(ss(s):se(s),:)=g;ll=ll+L;end
end

function state=local_viterbi(Y,ss,se,m)
logB=local_logB(Y,m.mu,m.var);K=size(m.A,1);state=zeros(size(Y,1),1);lp=log(max(m.pi,realmin));lA=log(max(m.A,realmin));
for s=1:numel(ss),a=ss(s);b=se(s);T=b-a+1;d=zeros(T,K);psi=zeros(T,K);d(1,:)=lp'+logB(a,:);
 for t=2:T,for k=1:K,[v,j]=max(d(t-1,:)+lA(:,k)');d(t,k)=v+logB(a+t-1,k);psi(t,k)=j;end,end
 [~,q]=max(d(T,:));path=zeros(T,1);path(T)=q;for t=T-1:-1:1,path(t)=psi(t+1,path(t+1));end;state(a:b)=path;
end
end

function D=local_load(cfg)
labels={'Dia','Noche','Madrugada'};if ~isscalar(cfg.regime)||~ismember(cfg.regime,1:3),error('Seleccion invalida.');end;reg=labels{cfg.regime};file=fullfile(cfg.root,sprintf('AGC4_Trigger_%s_69_Jornadas.xlsx',reg));if ~isfile(file),error('No se encontro: %s',file);end
sheets=sheetnames(file);if numel(sheets)~=69,warning('Se esperaban 69 hojas y se encontraron %d.',numel(sheets));end
inputs={'weather/air_temperature.outside','weather/relative_humidity.outside','weather/radiation_global','weather/wind_speed','weather/rain_state','weather/par.outside','weather/heat_emission','compartment/heating_lower_circuit/pipe_temperature','compartment/screen_energy/screen_position','compartment/screen_blackout/screen_position','compartment/water_supply/water_flow_duration','compartment/window_position_lee_side','compartment/window_position_wind_side','compartment/co2_actuation_state','compartment/lamps_activation_percentage'};
outputs={'compartment/water_drain/water_volume','compartment/water_drain/ec','compartment/air_temperature','compartment/relative_humidity','compartment/par','compartment/co2_concentration','compartment/substrate/relative_permittivity','compartment/substrate/bulk_ec','compartment/substrate/substrate_temperature','compartment/substrate/relative_permittivity.1','compartment/substrate/bulk_ec.1','compartment/substrate/substrate_temperature.1','compartment/substrate/relative_permittivity.2','compartment/substrate/bulk_ec.2','compartment/substrate/substrate_temperature.2','compartment/leaf_temperature','compartment/leaf_temperature.1'};
X=[];dayID=[];rowID=[];timeText={};
for s=1:numel(sheets),T=readtable(file,'Sheet',sheets{s},'VariableNamingRule','preserve');need=[{'hora'},inputs,outputs];miss=setdiff(need,T.Properties.VariableNames);if ~isempty(miss),error('Hoja %s: faltan columnas: %s',sheets{s},strjoin(miss,', '));end;M=NaN(height(T),numel(inputs)+numel(outputs));for j=1:numel(inputs),M(:,j)=double(T.(inputs{j}));end;for j=1:numel(outputs),M(:,numel(inputs)+j)=double(T.(outputs{j}));end;ir=find(strcmp(inputs,'compartment/water_supply/water_flow_duration'));d=[NaN;diff(M(:,ir))];d(d<0)=NaN;M(:,ir)=d;X=[X;M];dayID=[dayID;repmat(s,height(T),1)];rowID=[rowID;(1:height(T))'];h=string(T.hora);timeText=[timeText;cellstr(h)];end
inputs{strcmp(inputs,'compartment/water_supply/water_flow_duration')}='IrrigationIncrement';D=struct('file',file,'regime',reg,'X',X,'inputs',{inputs},'outputs',{outputs},'dayID',dayID,'rowID',rowID,'timeText',{timeText});
end

function [Z,med,scale,method,keep]=local_scale(X,names)
keep=true(1,size(X,2));med=NaN(1,size(X,2));scale=NaN(1,size(X,2));method=strings(1,size(X,2));for j=1:size(X,2),x=X(:,j);med(j)=median(x);q=prctile(x,[25 75]);iq=q(2)-q(1);sd=std(x);if isfinite(iq)&&iq>0,scale(j)=iq;method(j)="IQR";elseif isfinite(sd)&&sd>0,scale(j)=sd;method(j)="STD_FALLBACK";else,keep(j)=false;method(j)="CONSTANT_SKIPPED";warning('Variable constante omitida: %s',names{j});end,end;Z=(X(:,keep)-med(keep))./scale(keep);
end
function s=local_safe(s),s=regexprep(s,'[^A-Za-z0-9]+','_');s=regexprep(s,'^_+|_+$','');if strlength(s)>60,s=extractBefore(s,61);end,end
