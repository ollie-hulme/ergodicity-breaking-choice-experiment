function plotHLM(realData,runModelNum,whichPlots,flagPX,whichDir)
%% plotHLM
% plotHLM plots and post-processes HLM results, generating histograms, 
% computes MAPs, and generally visualises data
% It takes the following inputs:

% realData is a number indicating which data to plot.
% 1 is real agents
% 2 is synthetic agents for model recovery
% 3 is synthetic agetns for parameter recovery

% runModelNum
% 1 is parameter estimation via the condition specific isoelastic model (not strictly a latent mixture model)
% 2 is model comparison via a full HLM
% 3 is parameter estimation via a full HLM model

% whichPlots
% specifies which plots to generate according to the numbers listed in the
% section headings of this function. This should be a vector or scalar

% flagPX
% flags parameter expansion of the model indicator variable
% 1 indicates expansion upto a value of 12
% 0 indicates no expansion

%% Housekeeping
close all
cd '/mnt/projects/LogUtil/code/Behav/BayesianModeling/sandbox_ollie'

%% Set directories, figure properties and initial variables

figDir='/mnt/projects/LogUtil/code/Behav/BayesianModeling/sandbox_ollie/figs/';
if exist('realData','var')==0
    realData=input('real (n=18)(1) or synth data model recovery (n=27)(2) synth data param recovery (n=9) (3)');
end
if exist('runModelNum','var')==0
    runModelNum=input('which model number: 1 param estimation via iso; 2 model comparison; 3 param estimation via full model');
end
if exist('whichPlots','var')==0
    whichPlots=input('which plots to generate?');
end
if exist('whichDir','var')==0
    cd('samples_stats');nSamples=1;
else
    cd(['samples_stats/',whichDir]);sampList=dir('JAGS*');nSamples=numel(sampList);
end
switch realData
    case 1%real agents
        nSubjects=18;subjList=[1:4,6:19];%filters out subject 5, who is multiple SDs away
    case 2%synth agents for model recovery
        nSubjects=27; subjList=1:27;
    case 3%synth iso agents for param recovery across condition specific eta space
        nSubjects=9; subjList=1:9;
end

%% Set model-specific variables
switch runModelNum
    case 1
        modelName = 'Parameter estimation iso model';nModels=1;
    case 2
        modelName = 'Model selection full model';nModels=3;
    case 3
        modelName = 'Parameter estimation full model';nModels=3;
end

%% Figure properties
%set(0,'DefaultFigureWindowStyle','docked');
%eps = 0.01; binsc = eps/2:eps:1-eps/2; binse = 0:eps:1;

%% Loop over multiple JAGS samples results 
for samplp=1:nSamples
    
    if exist('whichDir','var')==0
    disp('load the jags output');uiopen;
    else
    load([sampList(samplp).name])   
    end

    %% 0 Plot Hyperpriors
    if find(whichPlots==0)
        
        figName='Hyperpriors';
        nFigs=nFigs+1;
        figure(nFigs);set(gcf,'color', 'w', 'units', 'normalized','paperpositionmode','auto');
        nRows=3;nCols=4;ct=1;
        
        subplot(nRows,nCols,ct);%Prior Log Mu Beta lin
        histogram(samples.lmu_beta_lin_m,nBins,'Normalization','probability');
        title('Prior LogMuBeta Lin', 'Interpreter', 'none');
        ct=ct+1;
        
        subplot(nRows,nCols,ct);%Target Prior Log Mu Beta lin
        lmu_beta_lin_unif = makedist('Uniform','lower',-2.3,'upper',1.61);
        lmu_beta_lin_unif_drawn = random(lmu_beta_lin_unif,size(samples.lmu_beta_lin));
        histogram(lmu_beta_lin_unif_drawn,nBins);
        title('Target Prior LogMuBeta Lin', 'Interpreter', 'none');
        ct=ct+1;
        
        subplot(nRows,nCols,ct);%Post Log Mu Beta lin add
        histogram(samples.lmu_beta_lin(:,:,1),nBins);
        title('Post LogMuBeta Lin Add', 'Interpreter', 'none');
        ct=ct+1;
        
        subplot(nRows,nCols,ct);%Post Log Mu Beta lin mult
        histogram(samples.lmu_beta_lin(:,:,2),nBins);
        title('Post LogMuBeta Lin Mult', 'Interpreter', 'none');
        ct=ct+1;
        
        subplot(nRows,nCols,ct);%Prior Log Sigma Beta lin
        histogram(samples.lsigma_beta_lin_m,nBins);
        title('Prior LogSigmaBeta Lin', 'Interpreter', 'none');
        ct=ct+1;
        
        subplot(nRows,nCols,ct);%Target Prior Log Sigma Beta lin
        lsigma_beta_lin_unif = makedist('Uniform','lower',0.01,'upper',1.13);
        lsigma_beta_lin_unif_drawn = random(lsigma_beta_lin_unif,size(samples.lsigma_beta_lin));
        histogram(lsigma_beta_lin_unif_drawn,nBins);
        title('Target Prior LogSigmaBeta lin', 'Interpreter', 'none');
        ct=ct+1;
        
        subplot(nRows,nCols,ct);%Post Log Sigma Beta lin
        histogram(samples.lsigma_beta_lin(:,:,1),nBins);
        title('Post LogSigmaBeta lin add', 'Interpreter', 'none');
        ct=ct+1;
        
        subplot(nRows,nCols,ct);%Post Log Sigma Beta lin
        histogram(samples.lsigma_beta_lin(:,:,2),nBins);
        title('Post LogSigmaBeta lin mult', 'Interpreter', 'none');
        ct=ct+1;
        
        subplot(nRows,nCols,ct);%Prior Log Tau Beta lin
        histogram(samples.ltau_beta_lin_m,nBins);xlim([0,200])
        title('Prior LogTauBeta lin', 'Interpreter', 'none');
        ct=ct+1;
        
        subplot(nRows,nCols,ct);%Target Prior Log Tau Beta lin
        ltau_beta_lin_unif_drawn = lsigma_beta_lin_unif_drawn.^-2;
        histogram(ltau_beta_lin_unif_drawn,nBins);xlim([0,200])
        title('Target Prior LogTauBeta lin', 'Interpreter', 'none');
        ct=ct+1;
        
        subplot(nRows,nCols,ct);%Post Log Tau Beta lin
        histogram(samples.ltau_beta_lin(:,:,1),nBins);xlim([0,200])
        title('Post LogTauBeta lin add', 'Interpreter', 'none');
        ct=ct+1;
        
        subplot(nRows,nCols,ct);%Post Log Tau Beta lin
        histogram(samples.ltau_beta_lin(:,:,2),nBins);xlim([0,200])
        title('Post LogTauBeta lin mult', 'Interpreter', 'none');
        
        figtitle(figName)
        saveas(gcf,[figDir,modelName,'_',figName,datestr(clock)],'epsc')
    end
    
    %% 1 Plot distributions of every variable in samples
    if find(whichPlots==1)
        figure
        sfields=fields(samples);nfields=length(sfields);
        for lp=1:nfields
            plotdim=ceil(sqrt(nfields));
            subplot(plotdim,plotdim,lp);
            eval(['histogram(','samples.',sfields{lp},'(:),100)']);
            ntitle(sfields{lp},'Interpreter', 'none');
        end
    end
    
    %% 2 Plot trace of every variable in samples
    if find(whichPlots==2)
        for lp=1:nfields
            % sort out dimensions
            eval(['x=size(','samples.',sfields{lp},');']);
            ndims=length(x);nchains=x(1);nsamps=x(2);try nsubs=x(3);end;try nconds=x(4);end
            
            %plot according to dimensions in the field
            switch ndims
                case 2
                    figure,eval(['plot(samples.',sfields{lp},'(:,:)'')']);figtitle(sfields{lp},'Interpreter', 'none');
                    figure,eval(['plot(log(samples.',sfields{lp},'(:,:)''))']);figtitle(['log_',sfields{lp}],'Interpreter', 'none');
                case 3
                    nsubs=x(3);plotdim=ceil(sqrt(nsubs));
                    figure
                    for i=1:nsubs
                        subplot(plotdim,plotdim,i);eval(['plot(samples.',sfields{lp},'(:,:,',num2str(i),')'')']);ntitle(num2str(i),'Interpreter', 'none');
                    end
                    figtitle(sfields{lp},'Interpreter', 'none');
                    figure
                    for i=1:nsubs
                        subplot(plotdim,plotdim,i);eval(['plot(log(samples.',sfields{lp},'(:,:,',num2str(i),'))'')']);ntitle(num2str(i),'Interpreter', 'none');
                    end
                    figtitle(['log_',sfields{lp}],'Interpreter', 'none');
                case 4
                    nsubs=x(3);nconds=x(4);plotdim=ceil(sqrt(nsubs));
                    for j=1:nconds
                        figure
                        for i=1:nsubs
                            subplot(plotdim,plotdim,i)
                            eval(['plot(samples.',sfields{lp},'(:,:,',num2str(i),',',num2str(j),')'')']);
                        end
                        figtitle([sfields{lp},'cond',num2str(j)],'Interpreter', 'none');
                        
                        figure
                        for i=1:nsubs
                            subplot(plotdim,plotdim,i)
                            eval(['plot(log(samples.',sfields{lp},'(:,:,',num2str(i),',',num2str(j),'))'')']);
                        end
                        figtitle(['log_',sfields{lp},'cond',num2str(j)],'Interpreter', 'none');
                    end
            end
        end
    end
    
    %% 3 Plot choice proportions
    if find(whichPlots==3)
        
        %% load data 
        load data_ergChoiceProportions
        
        %% rain clouds for proportions 
        figName='choicePropotions raincloud';figure; boxHeight = 0.4;transp=0.3;widthLine=0.5; %general constants      
        hold on;ylim([-1.75,3.5]);xlim([0,1]); yRainPlots=[-0.5,-1.25];colors=['b','r','k'];yticks(0:3);xticks(0:0.5:1);axis square  
        for lp=1:2           
        X{lp}=ergChoiceAllData{:,lp};yRainPlot=yRainPlots(lp);%mean level of rain plot and set data and color of data   
        [a,b] = ksdensity(X{lp},'Support',[0,1]); %plot density      
        h1 = area(b,a); set(h1, 'FaceColor', colors(lp),'FaceAlpha',transp,'EdgeAlpha',0);      
        if lp==1
        jit = (rand(size(X{lp})) - 0.5) * boxHeight;% jitter for raindrops
        end
        quants = quantile(X{lp},[0.25 0.75 0.5 0.02 0.98]);% info for making boxplot    
        x1=quants(1);x2=quants(2);y1=yRainPlot-(boxHeight*0.5);y2=yRainPlot+(boxHeight*0.5);
        x=[x1,x1,x2,x2]; y=[y1,y2,y2,y1];%set coordinates of patch
        h2=patch(x,y,'r');set(h2,'EdgeColor','k','LineWidth',widthLine,'FaceColor',colors(lp),'FaceAlpha',transp,'EdgeAlpha',0);
        h3 = line([quants(3) quants(3)],[yRainPlot+0.5*boxHeight yRainPlot-0.5*boxHeight],'col','k','LineWidth',widthLine);% mean line
        h4 = line([quants(2) quants(5)],[yRainPlot yRainPlot],'col','k','LineWidth',widthLine);h5 = line([quants(1) quants(4)],[yRainPlot yRainPlot],'col','k','LineWidth',widthLine);% whiskers
        h5 = scatter(X{lp},jit + yRainPlot,'SizeData',20,'MarkerFaceColor',colors(lp),'MarkerEdgeColor','none','MarkerFaceAlpha',transp,'MarkerEdgeAlpha',transp);% raindrops
        Y{lp}=jit + yRainPlot;
        end
        h6=patchline([0.5,0.5],[-4,3.5],'EdgeColor','k','LineStyle','--','EdgeAlpha',transp); 
        saveas(gcf,[figDir,'_',figName,datestr(clock)],'epsc');

        %% rainclouds for delta choice proportions
        figure, figName='delta_choicePropotions raincloud';figure; 
        hold on;ylim([-1,1.5]);xlim([-1,1]); yRainPlots=[-0.5,-1.5];colors=['k'];yticks(0:2);xticks(-1:0.5:1);axis square      
        for lp=1
        X{lp}=ergChoiceAllData{:,2}-ergChoiceAllData{:,1};yRainPlot=yRainPlots(lp);%mean level of rain plot and set data and color of data            
        [a,b] = ksdensity(X{lp},'Support',[-1,1]); %plot density      
        h1 = area(b,a); set(h1, 'FaceColor', colors(lp),'FaceAlpha',transp,'EdgeAlpha',0);      
        jit = (rand(size(X{lp})) - 0.5) * boxHeight;% jitter for raindrops
        quants = quantile(X{lp},[0.25 0.75 0.5 0.02 0.98]);% info for making boxplot    
        x1=quants(1);x2=quants(2);y1=yRainPlot-(boxHeight*0.5);y2=yRainPlot+(boxHeight*0.5);
        x=[x1,x1,x2,x2]; y=[y1,y2,y2,y1];%set coordinates of patch
        h2=patch(x,y,'r');set(h2,'EdgeColor','k','LineWidth',widthLine,'FaceColor',colors(lp),'FaceAlpha',transp,'EdgeAlpha',0);
        h3 = line([quants(3) quants(3)],[yRainPlot+0.5*boxHeight yRainPlot-0.5*boxHeight],'col','k','LineWidth',widthLine);% mean line
        h4 = line([quants(2) quants(5)],[yRainPlot yRainPlot],'col','k','LineWidth',widthLine);h5 = line([quants(1) quants(4)],[yRainPlot yRainPlot],'col','k','LineWidth',widthLine);% whiskers
        h5 = scatter(X{lp},jit + yRainPlot,'SizeData',20,'MarkerFaceColor',colors(lp),'MarkerEdgeColor','none','MarkerFaceAlpha',transp,'MarkerEdgeAlpha',transp);% raindrops
        h6=patchline([0,0],[-4,3.5],'EdgeColor','k','LineStyle','--','EdgeAlpha',transp);
        end
        saveas(gcf,[figDir,'_',figName,datestr(clock)],'epsc');
       
        %% raincloud plots for euclidean distances from map estimates to model predictions
        load data_euclidDistancesEta
        figName='euclidean distances of eta to predictions raincloud';figure; boxHeight = 0.4;transp=0.3;widthLine=0.5; %general constants      
        hold on;ylim([-1.75,1.5]);xlim([0,3]); yRainPlots=[-0.5,-1.25];colors=['b','r','k'];yticks(0:3);xticks(0:5);axis square  
        for lp=1:2           
        X{lp}=euclidDistancesEta{:,lp};yRainPlot=yRainPlots(lp);%mean level of rain plot and set data and color of data   
        [a,b] = ksdensity(X{lp},'Support',[0,3]); %plot density      
        h1 = area(b,a); set(h1, 'FaceColor', colors(lp),'FaceAlpha',transp,'EdgeAlpha',0);      
        if lp==1
        jit = (rand(size(X{lp})) - 0.5) * boxHeight;% jitter for raindrops
        end
        quants = quantile(X{lp},[0.25 0.75 0.5 0.02 0.98]);% info for making boxplot    
        x1=quants(1);x2=quants(2);y1=yRainPlot-(boxHeight*0.5);y2=yRainPlot+(boxHeight*0.5);
        x=[x1,x1,x2,x2]; y=[y1,y2,y2,y1];%set coordinates of patch
        h2=patch(x,y,'r');set(h2,'EdgeColor','k','LineWidth',widthLine,'FaceColor',colors(lp),'FaceAlpha',transp,'EdgeAlpha',0);
        h3 = line([quants(3) quants(3)],[yRainPlot+0.5*boxHeight yRainPlot-0.5*boxHeight],'col','k','LineWidth',widthLine);% mean line
        h4 = line([quants(2) quants(5)],[yRainPlot yRainPlot],'col','k','LineWidth',widthLine);h5 = line([quants(1) quants(4)],[yRainPlot yRainPlot],'col','k','LineWidth',widthLine);% whiskers
        h5 = scatter(X{lp},jit + yRainPlot,'SizeData',20,'MarkerFaceColor',colors(lp),'MarkerEdgeColor','none','MarkerFaceAlpha',transp,'MarkerEdgeAlpha',transp);% raindrops
        Y{lp}=jit + yRainPlot;
        end
        for slp=1:18
        h6=patchline([X{1}(slp),X{2}(slp)],[Y{1}(slp),Y{2}(slp)],'EdgeColor','k','LineStyle','--','EdgeAlpha',transp); 
        end
        saveas(gcf,[figDir,'_',figName,datestr(clock)],'epsc');
        
        %% raincloud plots for differences in euclidean distances from map estimates to model predictions
        figName='delta euclidean distances of eta to predictions raincloud';figure; boxHeight = 0.4;transp=0.3;widthLine=0.5; %general constants      
        hold on;ylim([-1.75,5]);xlim([0,3]); yRainPlots=[-0.5,-1.25];colors=['k'];yticks(0:3);xticks(0:0.5:1);axis square  
        for lp=1         
        X{lp}=euclidDistancesEta{:,2}-euclidDistancesEta{:,1};yRainPlot=yRainPlots(lp);%mean level of rain plot and set data and color of data   
        [a,b] = ksdensity(X{lp}); %plot density      
        h1 = area(b,a); set(h1, 'FaceColor', colors(lp),'FaceAlpha',transp,'EdgeAlpha',0);      
        if lp==1
        jit = (rand(size(X{lp})) - 0.5) * boxHeight;% jitter for raindrops
        end
        quants = quantile(X{lp},[0.25 0.75 0.5 0.02 0.98]);% info for making boxplot    
        x1=quants(1);x2=quants(2);y1=yRainPlot-(boxHeight*0.5);y2=yRainPlot+(boxHeight*0.5);
        x=[x1,x1,x2,x2]; y=[y1,y2,y2,y1];%set coordinates of patch
        h2=patch(x,y,'r');set(h2,'EdgeColor','k','LineWidth',widthLine,'FaceColor',colors(lp),'FaceAlpha',transp,'EdgeAlpha',0);
        h3 = line([quants(3) quants(3)],[yRainPlot+0.5*boxHeight yRainPlot-0.5*boxHeight],'col','k','LineWidth',widthLine);% mean line
        h4 = line([quants(2) quants(5)],[yRainPlot yRainPlot],'col','k','LineWidth',widthLine);h5 = line([quants(1) quants(4)],[yRainPlot yRainPlot],'col','k','LineWidth',widthLine);% whiskers
        h5 = scatter(X{lp},jit + yRainPlot,'SizeData',20,'MarkerFaceColor',colors(lp),'MarkerEdgeColor','none','MarkerFaceAlpha',transp,'MarkerEdgeAlpha',transp);% raindrops
        Y{lp}=jit + yRainPlot;
        end
        %h6=patchline([0.5,0.5],[-4,3.5],'EdgeColor','k','LineStyle','--','EdgeAlpha',transp); 
        saveas(gcf,[figDir,'_',figName,datestr(clock)],'epsc');
        
        
        
        % histogram for proportions
        figName='choicePropotions frequency histogram';figure;
        histogram(ergChoiceAllData{:,1},0:.1:1,'FaceColor','b','FaceAlpha',.5,'EdgeColor','b');hold on
        histogram(ergChoiceAllData{:,2},0:0.1:1,'FaceColor','r','FaceAlpha',.5,'EdgeColor','r')
        axis square;saveas(gcf,[figDir,'_',figName,datestr(clock)],'epsc');
        
        
        
        
        figName='delta choicePropotions frequency histogram';
        figure
        deltaCP=ergChoiceAllData{:,2}-ergChoiceAllData{:,1};
        histogram(deltaCP,-1:.2:1,'FaceColor',[0.5,0.5,0.5],'EdgeColor','k')
        ylim([0,6]);
        axis square
        saveas(gcf,[figDir,'_',figName,datestr(clock)],'epsc');
        
        figName='delta choicePropotions raincloud';
        figure
        deltaCP=ergChoiceAllData{:,2}-ergChoiceAllData{:,1};
        raincloud_plots_raincloud_plot(deltaCP,[0.5,0.5,0.5])
        ylim([0,6]);
        axis square
        saveas(gcf,[figDir,'_',figName,datestr(clock)],'epsc');
        
    end
    
    %% Model specific plots
    switch runModelNum
        
        case 1 %Iso model conditionwise and subjectwise for parameter estimation
            
            %% Compile posterior distributions for group for mu_eta, sigma_eta
            nVals=25;% Set number of values to evaluate histograms over
            muEtaAdd=samples.mu_eta(:,:,1);muEtaAdd=muEtaAdd(:);%aggregate parameter over subjects
            muEtaMult=samples.mu_eta(:,:,2);muEtaMult=muEtaMult(:);
            sigEtaAdd=samples.sigma_eta(:,:,1);sigEtaAdd=sigEtaAdd(:);
            sigEtaMult=samples.sigma_eta(:,:,2);sigEtaMult=sigEtaMult(:);
            
            %% Compute frequencies at group level
            [freqsEta,valsEta]=hist3([muEtaAdd,muEtaMult],{linspace(-2,2,nVals) linspace(-2,2,nVals)});%freqs are absolute counts, vals are values of the param
            [freqsSig,valsSig]=hist3([sigEtaAdd,sigEtaMult],{linspace(0,5,nVals) linspace(0,5,nVals)});
            
            %% Compute MAPs at group level
            [pmapMuEta,indEta]=max(freqsEta(:));[x,y]=ind2sub(size(freqsEta),indEta);
            mapEtaVals=[valsEta{1}(x),valsEta{2}(y)];
            [pmapSigEta,indSig]=max(freqsSig(:));[x,y]=ind2sub(size(freqsSig),indSig);
            mapSigVals=[valsSig{1}(x),valsSig{2}(y)];
            
            %% 4 Plot joint distribution of eta mu for group
            if find(whichPlots==4)
                figure,hold on,figName='joint distribution of mu eta group';sanePColor(valsEta{1},valsEta{2},freqsEta'/sum(sum(freqsEta')));
                xlim([-2,2]),ylim([-2,2]),axis('square')
                line([0, 0], ylim, 'LineWidth', 2, 'Color', 'b','LineStyle','-');line(ylim,[1, 1], 'LineWidth', 1, 'Color', 'r','LineStyle','-')
                scatter(linspace(-2,2,nVals),linspace(-2,2,nVals),'k.');hold on
                colorbar
                figtitle(figName)
                saveas(gcf,[figDir,modelName,'_',figName,datestr(clock)],'epsc');
            end
            
            %% 5 Plot joint distribution of eta and beta for group
            % to be repeated once mu_beta is monitored
            if find(whichPlots==5)
            end
            
            %% 6 Plot joint distribution of sigma for group
            % to be replotted with more expansive prior allowing sigma to hit 5
            if find(whichPlots==6)
                figure,hold on,figName='joint distribution of sigma eta group';sanePColor(valsSig{1},valsSig{2},freqsSig'/sum(sum(freqsSig')));
                xlim([0,5]),ylim([0,5])
                scatter(linspace(0,5,nVals),linspace(0,5,nVals),'k.'),axis('square')
                colorbar
                figtitle(figName)
                saveas(gcf,[figDir,modelName,'_',figName,datestr(clock)],'epsc');
            end
            
            %% 7 Plot posterior eta from mu and sigma for group
            % this involves summing over the sample hyperpriors for mu & sigma of eta
            if find(whichPlots==7)
                nVals=11;
                figure,hold on,figName='posterior eta based computed from group posterior parameters';
                postEtaGroup=zeros(nVals);
                etaRange=linspace(-2,2,nVals);
                [X1,X2] = meshgrid(etaRange',etaRange');X = [X1(:) X2(:)];
                postEtaGroup=zeros(nVals);
                for i=1:length(muEtaAdd)%sum over pdfs for every pair of mu and sigma
                    postEtaGroup=postEtaGroup+reshape(mvnpdf(X,[muEtaAdd(i),muEtaMult(i)],[sigEtaAdd(i),sigEtaAdd(i)]),nVals,nVals);
                end
                postEtaGroup=postEtaGroup./(sum(sum(postEtaGroup)));
                sanePColor(etaRange,etaRange,postEtaGroup),axis('square'),hold on
                scatter(linspace(-2,2,nVals),linspace(-2,2,nVals),'k.');
                xlim([-2,2]);ylim([-2,2]);colorbar
                line([0, 0], ylim, 'LineWidth', 1, 'Color', 'b','LineStyle','-');line(ylim,[1, 1], 'LineWidth', 1, 'Color', 'r','LineStyle','-')
                figtitle(figName)
                saveas(gcf,[figDir,modelName,'_',figName,datestr(clock)],'epsc');
                
                %% Compute bayes factors for group
                postTimeStrong=postEtaGroup(ceil(0.75*nVals),ceil(0.5*nVals));%coordinates for the eta=0,1
                priorTimeStrong=1/(nVals*nVals);
                bf_time_null=postTimeStrong/priorTimeStrong;
                bf_time_iso_any=postTimeStrong./(fliplr(eye(nVals)).*postEtaGroup);%pairwise comparison with any on diagonal
                bf_time_iso_all=postTimeStrong./(sum(sum(fliplr(eye(nVals)).*postEtaGroup))/nVals);
                
            end
            
            %% 8 Plot posterior eta subjectwise concatenated
            if find(whichPlots==8)
                figure,figName='Posterior eta subject- & condition-wise: all subjects concatenated';hold on
                etaAdd=samples.eta(:,:,1:length(subjList),1);etaAdd=etaAdd(:);etaMult=samples.eta(:,:,1:length(subjList),2);etaMult=etaMult(:);%configure data
                hpdAdd=hpdi(etaAdd,95);hpdMult=hpdi(etaMult,95);%estimate highest posterior density cred intervals
                [f1,x1] = ksdensity(etaAdd,'NumPoints',1000);[f2,x2] = ksdensity(etaMult,'NumPoints',1000);%kernal density
                edges=-2.55:0.1:2.55;%set edges of histogram for eta
                histogram(etaAdd,edges,'Normalization','pdf');histogram(etaMult,edges,'Normalization','pdf')%generate histograms for each cond
                plot(x1,f1,'b','LineWidth', 1.5);plot(x2,f2,'r','LineWidth', 1.5);%plot kernal density estimates
                line([0, 0], ylim, 'LineWidth', 1.5, 'Color', 'b','LineStyle','--');line([1, 1], ylim, 'LineWidth', 1.5, 'Color', 'r','LineStyle','--');%plot fiducials
                mapAddMean=[x1(find(f1==max(f1)))];mapMultMean=[x2(find(f2==max(f2)))];%compute MAPs
                figtitle(figName);
                saveas(gcf,[figDir,modelName,'_',figName,datestr(clock)],'epsc');
            end
            
            %% 9 Posterior eta mu over group
            if find(whichPlots==9)
                figure,figName='Posterior eta mu subject- & condition-wise: all subjects concatenated';hold on
                etaAdd=samples.mu_eta(:,:,1);etaAdd=etaAdd(:);etaMult=samples.mu_eta(:,:,2);etaMult=etaMult(:);%configure data
                hpdAdd=hpdi(etaAdd,95);hpdMult=hpdi(etaMult,95);%estimate highest posterior density cred intervals
                [f1,x1] = ksdensity(etaAdd,'NumPoints',1000);[f2,x2] = ksdensity(etaMult,'NumPoints',1000);%kernal density
                edges=-2.55:0.1:2.55;%set edges of histogram for eta
                histogram(etaAdd,edges,'Normalization','pdf');histogram(etaMult,edges,'Normalization','pdf')%generate histograms for each cond
                plot(x1,f1,'b','LineWidth', 1.5);plot(x2,f2,'r','LineWidth', 1.5);%plot kernal density estimates
                line([0, 0], ylim, 'LineWidth', 1.5, 'Color', 'b','LineStyle','--');line([1, 1], ylim, 'LineWidth', 1.5, 'Color', 'r','LineStyle','--');%plot fiducials
                mapAddMean=[x1(find(f1==max(f1)))];mapMultMean=[x2(find(f2==max(f2)))];%compute MAPs
                figtitle(figName);
                saveas(gcf,[figDir,modelName,'_',figName,datestr(clock)],'epsc');
            end
            
            %% 10 Plot Posterior eta marginalised over subjects
            if find(whichPlots==10)
                nVals=100;
                figure, figName='Conditionwise posterior density eta marginalised over subjects';
                [freqsEta,valsEta]=hist3([etaAdd,etaMult],{linspace(-2,2,nVals),linspace(-2,2,nVals)});
                sanePColor(valsEta{1},valsEta{2},freqsEta'/sum(sum(freqsEta')))
                xlim([-2,2]),ylim([-2,2]),axis('square');hold on
                scatter(linspace(-2,2,nVals),linspace(-2,2,nVals),'k.')
                line([0, 0], ylim, 'LineWidth', 1, 'Color', 'b','LineStyle','-');line(ylim,[1, 1], 'LineWidth', 1, 'Color', 'r','LineStyle','-')
                if realData == 3
                    %[x,y]=meshgrid(-0.5:0.5:1.5,-.5:0.5:1.5);
                    hold on
                    %scatter(x(:),y(:),'r');
                    scatter([0,0.5,1,0,0,0.5],[0,0.5,1,0.5,1,1],'g')
                end
                colorbar
                figtitle(figName)
                saveas(gcf,[figDir,modelName,'_',figName,datestr(clock)],'epsc');
            end
            
            %% 11 Plot hypothetical distributions over eta space - time
            if find(whichPlots==11)
                figure, figName='Predicted eta distributions - time';
                mu = [0 1]; Sigma = [0.3,0.3];%time average hypothesis
                [X1,X2] = meshgrid(linspace(-2,2,nVals)', linspace(-2,2,nVals)');
                X = [X1(:) X2(:)];
                p = mvnpdf(X, mu, Sigma);
                sanePColor(valsEta{1},valsEta{2},reshape(p,nVals,nVals));
                xlim([-2,2]),ylim([-2,2]),axis('square');hold on
                line([0, 0], ylim, 'LineWidth', 1, 'Color', 'b','LineStyle','--');line(ylim,[1, 1], 'LineWidth', 1, 'Color', 'r','LineStyle','--')
                scatter(linspace(-2,2,nVals),linspace(-2,2,nVals),'k.')
                saveas(gcf,[figDir,modelName,'_',figName,datestr(clock)],'epsc');
                
                figure, figName='Predicted eta distributions - dynamic invariant';
                mu = [0 0]; Sigma = [100 99.9; 99.9 100];%time average hypothesis
                [X1,X2] = meshgrid(linspace(-2,2,nVals)', linspace(-2,2,nVals)');
                X = [X1(:) X2(:)];
                p = mvnpdf(X, mu, Sigma);
                sanePColor(valsEta{1},valsEta{2},reshape(p,nVals,nVals));
                xlim([-2,2]),ylim([-2,2]),axis('square');hold on
                line([0, 0], ylim, 'LineWidth', 1, 'Color', 'b','LineStyle','--');line(ylim,[1, 1], 'LineWidth', 1, 'Color', 'r','LineStyle','--')
                scatter(linspace(-2,2,nVals),linspace(-2,2,nVals),'k.')
                saveas(gcf,[figDir,modelName,'_',figName,datestr(clock)],'epsc');
            end
            
            %% 12 Plot Posterior eta single subjects
            if find(whichPlots==12)
                f=figure; figName='Posterior eta subjectwise & conditionwise: single subjects';
                f.PaperOrientation='Landscape';
                mapAdd=[];mapMult=[];
                for lp=1:length(subjList)%excludes outlier 5
                    subplot(3,6,lp),hold on;
                    etaAdd=samples.eta(:,:,lp,1);etaAdd=etaAdd(:);etaMult=samples.eta(:,:,lp,2);etaMult=etaMult(:);
                    [f1,x1] = ksdensity(etaAdd,'NumPoints',1000); [f2,x2] = ksdensity(etaMult,'NumPoints',1000);
                    edges=-2.55:0.1:2.55;
                    histogram(samples.eta(:,:,lp,1),edges,'Normalization','pdf');histogram(samples.eta(:,:,lp,2),edges,'Normalization','pdf')
                    plot(x1,f1,'b','LineWidth', 1),plot(x2,f2,'r','LineWidth', 1)
                    line([0, 0], ylim, 'LineWidth', 1, 'Color', 'b','LineStyle','--'),line([-1, -1], ylim, 'LineWidth', 1, 'Color', 'g','LineStyle','--'),line([1, 1], ylim, 'LineWidth', 1, 'Color', 'r','LineStyle','--')
                    xlim([-2.5,2.5]);ylim([0,6]);axis('square')
                    if lp>1
                        xticks([]),yticks([]);
                    end
                    mapAdd=[mapAdd;x1(find(f1==max(f1)))];mapMult=[mapMult;x2(find(f2==max(f2)))];
                    if realData==3
                        trueEtas=[-1,-1; 0 -1;1 -1;-1 0;0 0;1 0;-1 1;0 1;1 1];%Mult is first column, add is second
                        trueEtas=trueEtas(:,[2,1]);%add is first column mult is second
                        ntitle(num2str(trueEtas(lp,:)))
                    end
                end
                
                figtitle(figName);saveas(gcf,[figDir,modelName,'_',figName,datestr(clock)],'pdf');
            end
            
            %% 13 Plot MAP estimates of Eta
            if find(whichPlots==13)
                figure;figName='MAP eta';hold on,set(gca,'Color','b')
                sanePColor(valsEta{1},valsEta{2},freqsEta')
                scatter(mapAdd,mapMult,'g.');
                scatter(linspace(-2.5,2.5,100),linspace(-2.5,2.5,100),'k.');
                scatter(mapAddMean,mapMultMean,'m*')
                meanMapAdd=0.373; meanMapMult=1.374;meanMapAddCI=[0.188,0.558];meanMapMultCI=[1.210,1.537];  %stats from jasp analysis
                scatter(meanMapAdd,meanMapMult,'m.')
                errorbar(meanMapAdd,meanMapMult,meanMapAddCI(1)-meanMapAdd,meanMapAddCI(2)-meanMapAdd,...
                    meanMapMultCI(1)-meanMapMult,meanMapMultCI(2)-meanMapMult,'-m')
                xlim([-2.5,2.5]);ylim([-2.5,2.5])
                line([0, 0], ylim, 'LineWidth', 1, 'Color', 'b','LineStyle','--');line(ylim,[1, 1], 'LineWidth', 1, 'Color', 'r','LineStyle','--')
                axis('square');xlabel('eta additive');ylabel('eta multiplicative')
                figtitle(figName)
                if realData == 3
                    [x,y]=meshgrid(-1:0.5:1,-1:0.5:1),hold on, scatter(x(:),y(:),'r');
                end
                saveas(gcf,[figDir,modelName,'_',figName,datestr(clock)],'epsc');
            end
            
            %% 14 Plot Effect of condition on movement of MAP eta-beta for single subjects
            if find(whichPlots==14)
                mapEtaVals=[];mapBetaVals=[];day1=nan(18,2);day2=day1;
                for i=subjList%excludes outlier 5
                    etaAdd=samples.eta(:,:,i,1);etaAdd=etaAdd(:);
                    betaAdd=samples.beta_crra(:,:,i,1);betaAdd=betaAdd(:);
                    etaMult=samples.eta(:,:,i,2);etaMult=etaMult(:);
                    betaMult=samples.beta_crra(:,:,i,2);betaMult=betaMult(:);
                    [freqsEta,valsEta]=hist3([etaAdd,etaMult],[nVals,nVals]);%freqs are absolute counts, vals are values of the param
                    [freqsBeta,valsBeta]=hist3([betaAdd,betaMult],[nVals,nVals]);
                    [pmapMuEta,indEta]=max(freqsEta(:));[x,y]=ind2sub(size(freqsEta),indEta);
                    mapEtaVals(i,:)=[valsEta{1}(x),valsEta{2}(y)];
                    [pmapBeta,indBeta]=max(freqsBeta(:));[x,y]=ind2sub(size(freqsBeta),indBeta);
                    mapBetaVals(i,:)=[valsBeta{1}(x),valsBeta{2}(y)];
                    if i<10
                        day1(i,:)=[mapEtaVals(i,2),log(mapBetaVals(i,2))];%on day1 it is multiplicative condition
                        day2(i,:)=[mapEtaVals(i,1),log(mapBetaVals(i,1))];%on day2 it is additive condition
                    else
                        day1(i,:)=[mapEtaVals(i,1),log(mapBetaVals(i,1))];%on day1 it is add condition
                        day2(i,:)=[mapEtaVals(i,2),log(mapBetaVals(i,2))];%on day2 it is mult condition
                    end
                end
                
                figure;figName='beta-eta trajectories';subInds=[1:4,6:19];
                subplot(3,1,1);hold on;
                x=[mapEtaVals(subInds,1),mapEtaVals(subInds,2)];y=[log(mapBetaVals(subInds,1)),log(mapBetaVals(subInds,2))];
                for lp=1:length(subInds)
                    plot1=plot(x(lp,:),y(lp,:));plot1.Color(4)=0.4;plot1.LineWidth=2;
                    scatter(x(lp,1),y(lp,1),30,plot1.Color,'o','filled');scatter(x(lp,2),y(lp,2),30,plot1.Color,'o')
                end
                xlim([-2.5,2.5]);ylim([-5,8]);axis('square')
                
                subplot(3,1,2);hold on;subInds=10:19;
                x=[mapEtaVals(subInds,1),mapEtaVals(subInds,2)];y=[log(mapBetaVals(subInds,1)),log(mapBetaVals(subInds,2))];
                for lp=1:length(subInds)
                    plot1=plot(x(lp,:),y(lp,:),'r');plot1.Color(4)=0.4;plot1.LineWidth=2;
                    scatter(x(lp,1),y(lp,1),30,'r','o','filled');scatter(x(lp,2),y(lp,2),30,plot1.Color,'o')
                end
                xlim([-2.5,2.5]);ylim([-5,8]);axis('square')
                
                subplot(3,1,3),hold on;subInds=[1:4,6:9];
                x=[mapEtaVals(subInds,1),mapEtaVals(subInds,2)];y=[log(mapBetaVals(subInds,1)),log(mapBetaVals(subInds,2))];
                for lp=1:length(subInds)
                    plot1=plot(x(lp,:),y(lp,:),'b');plot1.Color(4)=0.4;plot1.LineWidth=2;
                    scatter(x(lp,1),y(lp,1),30,'b','o');scatter(x(lp,2),y(lp,2),30,plot1.Color,'o','filled')
                end
                xlim([-2.5,2.5]);ylim([-5,8]);axis('square'),figtitle(figName);
                saveas(gcf,[figDir,modelName,'_',figName,datestr(clock)],'pdf');
            end
            
            %% 15 Plot Eta versus TimeAv
            if find(whichPlots==15)
                figure,hold on,figName='Eta map estimates versus time averages';load('data_TimeAveragesChosenGambles')
                yyaxis left
                scatter(mapAdd,TimeAveragesChosenGambles(:,1),'b.')
                line([0, 0], ylim, 'LineWidth', 1, 'Color', 'b','LineStyle','--')
                ylabel('Time average additive growth rate')
                yyaxis right
                scatter(mapMult,TimeAveragesChosenGambles(:,2),'r.')
                ylabel('Time average multiplicative growth rate');xlabel('eta map')
                line([1, 1],ylim, 'LineWidth', 1, 'Color', 'r','LineStyle','--')
                axis('square')
                figtitle(figName)
                saveas(gcf,[figDir,modelName,'_',figName,datestr(clock)],'epsc');
            end
            
        case 2 %Model selection full model
            
            %% 16 Plot posterior z subjectwise
            if find(whichPlots==16)
                figure,figName='Posterior z subjectwise';zSum=zeros(18,3);%for summing over z values
                if exist('flagPX','var')==0
                flagPX=input('param expansion on z? 1 for yes (1+4+7+10 etc) 0 for no (1,2,3)');
                end
                for i=1:length(subjList)
                    try
                        numtrials=size(samples.z,1)*size(samples.z,2);
                        z1=length(find(samples.z(:,:,i)==1))/numtrials; z2=length(find(samples.z(:,:,i)==2))/numtrials; z3=length(find(samples.z(:,:,i)==3))/numtrials;
                        z4=length(find(samples.z(:,:,i)==4))/numtrials; z5=length(find(samples.z(:,:,i)==5))/numtrials; z6=length(find(samples.z(:,:,i)==6))/numtrials;
                        z7=length(find(samples.z(:,:,i)==7))/numtrials; z8=length(find(samples.z(:,:,i)==8))/numtrials; z9=length(find(samples.z(:,:,i)==9))/numtrials;
                        z10=length(find(samples.z(:,:,i)==10))/numtrials; z11=length(find(samples.z(:,:,i)==11))/numtrials; z12=length(find(samples.z(:,:,i)==12))/numtrials;
                        if flagPX==0
                            stats.zprops{i}=[z1,z2,z3];
                        elseif flagPX==1
                            stats.zprops{i}=[z1+z4+z7+z10,z2+z5+z8+z11,z3+z6+z9+z12];
                        end
                        subplot(nSubjects,1,i),bar(stats.zprops{i});ylim([0,1]),xlim([0.5,3.5]),ylabel(num2str(subjList(i)))
                        set(gca,'xtick',[]);set(gca,'ytick',[])
                        zSum(i,:)=stats.zprops{i};
                    catch
                    end
                end
                figtitle(figName)
                saveas(gcf,[figDir,modelName,'_',figName,datestr(clock)],'epsc')
                figName='Posterior z subjectwise as image';
                figure,imagesc(zSum),colormap(flipud(bone));
            end
            
            %% 17 Plot posterior z summed over group
            if find(whichPlots==17)
                figure,figName='Posterior z subjectwise summed over group';
                sumProbs=sum(zSum,1)/nSubjects;
                bar(sumProbs),figtitle(figName),xlim([0.5,3.5]),ylim([0,1])
                saveas(gcf,[figDir,modelName,'_',figName,datestr(clock)],'epsc')
                bf_posteriorZ=sumProbs'.*(1./sumProbs)
            end
            
            %% 18 Plot protected exceedance probabilities
            if find(whichPlots==18)
                dims=size(samples.z);%get dimensions of samples
                nsamps=prod(dims(1:2));%find n samples per subject over chains
                zSum(zSum==0)=1/(nsamps+1);%set 0 to upper bound probabilities assuming 1 observation in the next sample
                [post,out]=VBA_groupBMC(log(zSum)');%calculate pxps
                
                figure,figName='Protected exceedance probabilities';
                bar(out.pxp),xlim([0.5,3.5]),ylim([0,1]);
                saveas(gcf,[figDir,modelName,'_',figName,datestr(clock)],'epsc')
                bf_pxp=out.pxp'.*(1./out.pxp)
                
                figName='Model attributions';
                figure,imagesc(post.r'),colormap(flipud(bone)),axis off
                saveas(gcf,[figDir,modelName,'_',figName,datestr(clock)],'epsc')
                
                figName='Model attributions marginalised over subjects';
                marginalAttributions=sum(post.r')/sum(sum(post.r'));
                bar(marginalAttributions),xlim([0.5,3.5]),ylim([0,1]);
                saveas(gcf,[figDir,modelName,'_',figName,datestr(clock)],'epsc')
                bf_marginalAttributions=marginalAttributions'.*(1./marginalAttributions)
            end
            
            %% 19 Plot traces conditioned on z 
            if find(whichPlots==19)      
                figure;               
                for mlp=0:2
                    subplot(1,3,1+mlp),hold on
                    for ch=1:4
                        x=samples.eta_iso(ch,:,:);
                        plot(x((samples.z(ch,:,:)==1+mlp | samples.z(ch,:,:)== 4+mlp | samples.z(ch,:,:)== 7+mlp |samples.z(ch,:,:)==10+mlp))),hold on
                    end
                end, ntitle('eta_iso') 
                figure;               
                for mlp=0:2
                    subplot(1,3,1+mlp),hold on
                    for ch=1:4
                        x=samples.eta_tw(ch,:,:);
                        plot(x((samples.z(ch,:,:)==1+mlp | samples.z(ch,:,:)== 4+mlp | samples.z(ch,:,:)== 7+mlp |samples.z(ch,:,:)==10+mlp))),hold on
                    end
                end, ntitle('eta_tw') 
                 figure;              
                for mlp=0:2
                    subplot(1,3,1+mlp),hold on
                    for ch=1:4
                        x=samples.alphaGain(ch,:,:);
                        plot(x((samples.z(ch,:,:)==1+mlp | samples.z(ch,:,:)== 4+mlp | samples.z(ch,:,:)== 7+mlp |samples.z(ch,:,:)==10+mlp))),hold on
                    end
                    end,ntitle('alphaGain') 
                figure;               
                for mlp=0:2
                    subplot(1,3,1+mlp),hold on
                    for ch=1:4
                        x=samples.alphaLoss(ch,:,:);
                        plot(x((samples.z(ch,:,:)==1+mlp | samples.z(ch,:,:)== 4+mlp | samples.z(ch,:,:)== 7+mlp |samples.z(ch,:,:)==10+mlp))),hold on
                    end
                end,ntitle('alphaLoss')  
                figure;              
                for mlp=0:2
                    subplot(1,3,1+mlp),hold on
                    for ch=1:4
                        x=samples.lambda(ch,:,:);
                        plot(x((samples.z(ch,:,:)==1+mlp | samples.z(ch,:,:)== 4+mlp | samples.z(ch,:,:)== 7+mlp |samples.z(ch,:,:)==10+mlp))),hold on
                    end
                end, ntitle('lambda')             
            end
            
        case 3 %Parameter estimation full model
            
            modInd=input('which model 1=time, 2=pt, 3=iso?');
            
            switch modInd
                case 1             
                case 2
                    figure
                    subplot(1,4,1),histogram(samples.lambda(:))
                    subplot(1,4,2),histogram(samples.alphaGain(:))
                    subplot(1,4,3),histogram(samples.alphaLoss(:))
                    subplot(1,4,4),histogram(samples.beta_pt(:,:,:,1),0:.05:1,'FaceColor','b','FaceAlpha',.5,'EdgeColor','b'),hold on
                    subplot(1,4,4),histogram(samples.beta_pt(:,:,:,2),0:.05:1,'FaceColor','r','FaceAlpha',.5,'EdgeColor','r'),xlim([0,1])
                case 3
            end
    end
end %end of sample loop (samplp)