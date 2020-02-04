function decoder = dffPxlDecoder(dff,outcome,method,cfg)
% decoder = dffPxlDecoder(rec,spockFlag,decodeWhat,method,correctViewAng,derivFlag,cfg,ROIflag,ROIlbl)
% linear decoder of task variables based on simulatenous activity of pixels
% or ROIs
%
% input:
%   dff: trials x ROIs x time
%   outcome: trials x 1 (populated with -1 or 1s)
%   method: 'ridge' (default) or 'lasso'
%   cfg: structure with analysis config, refer to bottom of function. If
%        empty will be filled with defaults (recommended)
%
% output is a data structure with analysis results
%
% 
% Original code found here: https://github.com/BrainCOGS/widefieldImaging

%% defaults, paths etc
if nargin < 2 || isempty(outcome); error('not enough input arguments'); end
if nargin < 3 || isempty(method);         method     = 'ridge';         end
if nargin < 4 || isempty(cfg);            cfg        = struct([]);      end


% cfg(1).decodeWhat       = decodeWhat;

switch method
  case 'ridge'
    cfg(1).method    = method;
  case 'lasso'
    cfg(1).method    = method;
    cfg(1).l2Penalty = 0;
  case 'elasticNet'
    cfg(1).method    = 'lasso';
    cfg(1).l2Penalty = 1;
end

cfg               = populateCfg(cfg);



%% set up response vector

y = outcome;

%% set up predictor matrix 
% [nTrials,nTime,nROI]   = size(dffTrials);
% 
% % if one decoder per time point, each slice of a 3d matrix is a predictor
% % matrix for that time point (trials x pixels)
% predMat              = nan(nTrials,nROI,nTime);
% for iTime = 1:nTime
%     predMat(:,:,iTime) = squeeze(dffTrials(:,iTime,:));
% end

[nTrials,nROI,nTime]   = size(dff);
predMat = dff;

% clean up nans
nanidx              = isnan(y(:,1));
y(nanidx,:)         = [];
predMat(nanidx,:,:) = [];

% compile info
decoder.nPoints     = nTime;
decoder.nTrials     = nTrials;
decoder.nPxls       = nROI;
decoder.cfg         = cfg;


%% set up xval


warning('off', 'stats:cvpartition:KFoldMissingGrp');
randGenerator         = RandStream('mt19937ar', 'Seed', 723176);
pseudoExp             = cvpartition(y(:,end), 'KFold', cfg.numCVFolds, randGenerator);
for iMC = 2:cfg.numCVSamples
  pseudoExp(iMC)      = repartition(pseudoExp(1));
end
warning('on', 'stats:cvpartition:KFoldMissingGrp');

%% train / test decoder

jitter                = cfg.stateJitter;
cl                    = cfg.confidLevel;
l2                    = cfg.l2Penalty;
lambda                = cfg.lambdaRatio;
method                = cfg.method;

bestCoeff             = cell(1,nTime);
accuracy              = cell(1,nTime);
accInterval           = cell(1,nTime);
coeff                 = cell(1,nTime);
fitInfo               = cell(1,nTime);

tic;

if ~cfg.singleDecoder
  parfor iTime = 1:nTime  % #ok<PFUIX> %PARFOR
    if size(y,2) > 1
      target          = y(:,iTime);
    else
      target          = y;
    end
    
    % add jitter to all-zero pixels 
    pred     = predMat(:,:,iTime);
    allzeros = sum(pred==0) == size(pred,1);
    pred(:,allzeros) = pred(:,allzeros) + jitter(1) * randn(size(pred(:,allzeros)));
    
    switch method
      case 'lasso'
        [coeff{iTime}, fitInfo{iTime}] ...
                       = elasticNetRegression( pred, target, [], ...
                                              {'LeastR',jitter}, l2, lambda, pseudoExp, [], [] ); 
        fitInfo{iTime} = rmfield(fitInfo{iTime},{'X','y','w'});
        
      case 'ridge'
        [coeff{iTime}, fitInfo{iTime}] ...
                       = ridgeRegression( pred, target, lambda, pseudoExp, [], jitter); 
    end
    
    bestCoeff{iTime}  = coeff{iTime}(:,fitInfo{iTime}.IndexMinMSE);
    
    % Compute accuracy as the fraction of correct predictions, with spreads across CV samples
    predAccuracy      = nan(size(fitInfo{iTime}.CVTestPrediction));

    % Cross-validated prediction accuracy
    for iMC = 1:size(fitInfo{iTime}.CVTestPrediction,1)
      for iFold = 1:size(fitInfo{iTime}.CVTestPrediction,2)
        % Compute optimal threshold for separating training set populations according to truth
        trainSel      = fitInfo{iTime}.CVExperiments(iMC).training(iFold);


          prediction  = fitInfo{iTime}.CVTrainPrediction{iMC,iFold}(:,fitInfo{iTime}.Index1SE);
          try
            threshold   = optimalClassThreshold(prediction, target(trainSel)>0);
          catch
            threshold   = getThreshold(prediction,0);
          end
          % Compute accuracy by applying threshold to prediction in test set and comparing to truth
          testSel     = fitInfo{iTime}.CVExperiments(iMC).test(iFold);
          prediction  = fitInfo{iTime}.CVTestPrediction{iMC,iFold}(:,fitInfo{iTime}.Index1SE);
          predTruth   = target(testSel);
          prediction(prediction < threshold) = -1;
          prediction(prediction > threshold) =  1;
          predAccuracy(iMC,iFold) = sum(prediction == predTruth)./numel(predTruth); 

      end
    end

    accuracy{iTime}    = nanmean(predAccuracy(:));
    accInterval{iTime} = quantile(predAccuracy(:), cl); 

  end
else
  error('single decoder not yet implemented')
end



%% Consolidate values
decoder.allCoeffs  = coeff;
decoder.fitInfo    = fitInfo;
decoder.accuracy   = cell2mat(accuracy);
decoder.weights    = cell2mat(bestCoeff)';

for iTime = 1:nTime
  decoder.accuracyCI(iTime,:) = accInterval{iTime};
end

% weights in brain image x time points (3d matrix)
% decoder.weights  = conditionDffMat(decoder.weights,bc,br,[nX nY nTime]);

clear bestCoeff accuracy coeff 

%% do shuffling if necessary
if cfg.numShuffles > 1
  

  bestCoeff             = cell(cfg.numShuffles,iTime);
  accuracy              = cell(cfg.numShuffles,iTime);
  numShuff              = cfg.numShuffles;
  numCVFolds            = cfg.numCVFolds;
  numCVSamples          = 5; % in each shuffle just do 5 runs of xval
  
  % use best lambda from actual data
  bestLambdas           = zeros(1,iTime);
  for iTime = 1:nTime
    bestLambdas(iTime)  = lambda([fitInfo{iTime}.IndexMinMSE]);
  end

  parfor iTime = 1:nTime
    if size(y,2) > 1
      target          = y(:,iTime);
    else
      target          = y;
    end
    
    % add jitter to all-zero pixels 
    pred     = predMat(:,:,iTime);
    allzeros = sum(pred==0) == size(pred,1);
    pred(:,allzeros) = pred(:,allzeros) + jitter(1) * randn(size(pred(:,allzeros))); %#ok<PFBNS>
    
    for iShuff = 1:numShuff 
      
      iTarget         = target(randperm(numel(target)));
      
      switch method 
        case 'lasso'
          [bestCoeff{iShuff,iTime}, sfitInfo] = elasticNetRegression( pred, iTarget, [],                         ...
                                                                      {'LeastR',jitter}, l2, bestLambdas(iTime), ...
                                                                      numCVFolds, numCVSamples, [] );
          
        case 'ridge'
          [bestCoeff{iShuff,iTime}, sfitInfo] = ridgeRegression( pred, iTarget, bestLambdas(iTime), ...
                                                                 numCVFolds, numCVSamples, jitter);
          
      end

      % Compute accuracy as the fraction of correct predictions, with spreads across CV samples
      predAccuracy      = nan(size(sfitInfo.CVTestPrediction));
      
      % Cross-validated prediction accuracy
      for iMC = 1:size(sfitInfo.CVTestPrediction,1)
        for iFold = 1:size(sfitInfo.CVTestPrediction,2)
          % Compute optimal threshold for separating training set populations according to truth
          trainSel      = sfitInfo.CVExperiments(iMC).training(iFold);
          

            prediction  = sfitInfo.CVTrainPrediction{iMC,iFold};
            try
              threshold   = optimalClassThreshold(prediction, target(trainSel)>0);
            catch
              threshold   = getThreshold(prediction,0);
            end
            % Compute accuracy by applying threshold to prediction in test set and comparing to truth
            testSel     = sfitInfo.CVExperiments(iMC).test(iFold);
            prediction  = sfitInfo.CVTestPrediction{iMC,iFold};
            predTruth   = target(testSel);
            prediction(prediction < threshold) = -1;
            prediction(prediction > threshold) =  1;
            predAccuracy(iMC,iFold) = sum(prediction == predTruth)./numel(predTruth); 


        end
      end
      
      accuracy{iShuff,iTime}    = nanmean(predAccuracy(:));
    end
  end
  
  decoder.shuffle.coeffs        = cell2mat(bestCoeff);
  decoder.shuffle.accuracy      = cell2mat(accuracy);
  

end

%% shut down parallel pool
% delete(poolobj);

%% save 
% save(fn,'decoder','cfg','-v7.3')

% %% plot decoding accuracy
% if ~isempty(dir([fn '.pdf'])); delete([fn '.pdf']); end
% 
% switch cfg.timeOrSpace
%   case 'time'
%     taxis = cfg.timeBins;
%     xlbl  = 'Time (s)';
%   case 'space'
%     taxis = cfg.posBins;
%     xlbl  = 'y position (cm)';
% end
% 
% if ~ROIflag
%   figure; wf.applyFigDefaults(gcf,[2 1],'w'); hold on
%   plot(taxis, decoder.accuracy, '-', 'color', 'k', 'linewidth', 1.5)
%   plot(taxis, decoder.accuracyCI(:,1), '--', 'color', 'k',   'linewidth', .75)
%   plot(taxis, decoder.accuracyCI(:,2), '--', 'color', 'k',   'linewidth', .75)
% 
%   if isfield(decoder,'shuffle')
%     m   = nanmean(decoder.shuffle.accuracy);
%     sem = nanstd(decoder.shuffle.accuracy)./sqrt(cfg.numShuffles);
% 
%     plot(taxis, m, '-', 'color', wf.lightgray, 'linewidth', 1)
%     plot(taxis, m+sem, '--', 'color', wf.lightgray,   'linewidth', .5)
%     plot(taxis, m-sem, '--', 'color', wf.lightgray,   'linewidth', .5)
%   end
% 
%   wf.applyAxisDefaults(gca,'k');
%   wf.applyAxisLbls(gca,xlbl,'Decoding Acc.',sprintf('%s - %s trials',cfg.decodeWhat,cfg.trialType))
%   saveas(gcf,[fn '_accuracy'])
%   export_fig([fn '.pdf'],'-append')
%   close
% 
%   %% plot decoding weights
%   [nr,nc]   = subplotOrg(nTime,8);
%   cmap      = colormap(red2blue);
%   cmap(1,:) = [0 0 0];
%   figure; wf.applyFigDefaults(gcf,[nc nr],'k'); hold on
%   for iTime = 1:nTime
%     subplot(nr,nc,iTime)
%     thisweight = decoder.weights(:,:,iTime);
%     maxc       = max(abs(thisweight(:)));
%     if isnan(maxc); continue; end
%     if maxc == 0; maxc = 0.01; end
%     imagesc(thisweight,[-maxc maxc]); colormap(cmap)
%     axis image; axis off
%     title(['y = ' num2str(taxis(iTime))],'color','w','fontsize',14,'fontweight','bold')
%   end
%   saveas(gcf,[fn '_weights'])
%   export_fig([fn '.pdf'],'-append')
%   close
% else
%   %%
%   figure; wf.applyFigDefaults(gcf,[4 2],'w'); 
%   
%   subplot(2,4,[1 2]); hold on
%   plot(taxis, decoder.accuracy, '-', 'color', 'k', 'linewidth', 1.5)
%   plot(taxis, decoder.accuracyCI(:,1), '--', 'color', 'k',   'linewidth', .75)
%   plot(taxis, decoder.accuracyCI(:,2), '--', 'color', 'k',   'linewidth', .75)
% 
%   if isfield(decoder,'shuffle')
%     m   = nanmean(decoder.shuffle.accuracy);
%     sem = nanstd(decoder.shuffle.accuracy)./sqrt(cfg.numShuffles);
% 
%     plot(taxis, m, '-', 'color', wf.lightgray, 'linewidth', 1)
%     plot(taxis, m+sem, '--', 'color', wf.lightgray,   'linewidth', .5)
%     plot(taxis, m-sem, '--', 'color', wf.lightgray,   'linewidth', .5)
%   end
% 
%   wf.applyAxisDefaults(gca,'k');
%   wf.applyAxisLbls(gca,xlbl,'Decoding Acc.',sprintf('%s - %s trials',cfg.decodeWhat,cfg.trialType))
%   
%   subplot(2,4,[5 6]); hold on
%   imagesc(taxis,1:nROI,decoder.weights',[-abs(max(decoder.weights(:))) abs(max(decoder.weights(:)))]); 
%   axis tight; colormap red2blue
%   colorbar('location','southoutside','position',[.1 .04 .1 .02])
%   set(gca,'ytick',1:nROI,'yticklabel',ROIlbl)
%   wf.applyAxisDefaults(gca,'k'); 
%   wf.applyAxisLbls(gca,xlbl,[],'Decoding weights (a.u.)')
%   
%   subplot(2,4,[3 4]); hold on
%   colors = jet(nROI);%feval(widefieldParams.colors,nROI);
%   weightmean = mean(decoder.weights);
%   for iROI = 1:nROI
%     bar(iROI,weightmean(iROI),'edgecolor',colors(iROI,:),'facecolor',colors(iROI,:)); 
%   end
%   set(gca,'xtick',1:nROI,'xticklabel',ROIlbl)
%   rotateXLabels(gca,30)
%   wf.applyAxisDefaults(gca,'k'); axis tight
%   wf.applyAxisLbls(gca,[],'Decoding weight (a.u.)','Average weight')
%   
%   subplot(2,4,[7 8]); hold on
%   weightmean = decoder.weights(decoder.accuracy == max(decoder.accuracy),:);
%   for iROI = 1:nROI
%     bar(iROI,weightmean(iROI),'edgecolor',colors(iROI,:),'facecolor',colors(iROI,:)); 
%   end
%   set(gca,'xtick',1:nROI,'xticklabel',ROIlbl)
%   rotateXLabels(gca,30)
%   wf.applyAxisDefaults(gca,'k'); axis tight
%   wf.applyAxisLbls(gca,[],'Decoding weight (a.u.)','Weight at max accuracy')
% end
% 
% catch ME
%   displayException(ME)
% end
end



%% --------  compile cfg
function cfg = populateCfg(cfg)

%% decoder config
if ~isfield(cfg,'trialType')
  cfg(1).trialType     = 'all'; % choice, prev. choice, evidence, position, view angle
end
% if ~isfield(cfg,'decodeWhat')
%   cfg(1).decodeWhat    = 'choice'; % choice, prev. choice, evidence, position, view angle
% end
if ~isfield(cfg,'singleDecoder')
  cfg(1).singleDecoder = false; % true to train one decoder for all time points
end
if ~isfield(cfg,'singleDecoderTimePoint')
  cfg(1).singleDecoderTimePoint = 300; % if single decoder, which point to train it on?
end
if ~isfield(cfg,'singleCut')
  cfg(1).singleCut     = false; % if false will use different categorization boundaries per time point
end

if ~isfield(cfg,'spatialBinFactor')
  cfg(1).spatialBinFactor = 4; % spatial binning
end


%% optimization config
if ~isfield(cfg,'method')
  cfg(1).method        = 'ridge'; % ridge or lasso (if the latter setting l2penalty to true will do elastic net)
end
if ~isfield(cfg,'lambdaRatio')
  switch cfg.method
    case 'lasso'
      cfg(1).lambdaRatio   = flip(10.^(linspace(-3, log10(0.5), 15))); % l1/l2 penalty
    case 'ridge'
      cfg(1).lambdaRatio   = flip(10.^(linspace(-3, log10(1), 20))); % l1/l2 penalty
  end
end
if ~isfield(cfg,'stateJitter')
  cfg.stateJitter      = [1e-6,1e-5]; % to prevent degeneracies in population states e.g. induced by zero cue-locked amplitudes
end
if ~isfield(cfg,'l2penalty')
  cfg(1).l2Penalty     = 0; % Elastic net vs. lasso
end
if ~isfield(cfg,'xvalNFold')
  cfg(1).numCVFolds    = 3;
end
if ~isfield(cfg,'xvalNSamples')
  cfg(1).numCVSamples  = 50;
end
if ~isfield(cfg,'numShuffles')
  cfg(1).numShuffles   = 10;
end
if ~isfield(cfg,'confidLevel')
  cfg(1).confidLevel   = normcdf([-1 1], 0, 1);
end

end

%%
function threshold = getThreshold(prediction,cut)

threshold = min(prediction) + (max(prediction)-min(prediction))/2;

if cut == 0
  if abs(threshold - cut) > .1;     threshold = cut; end
else
  if abs(threshold - cut) > .1*cut; threshold = cut; end
end

end

function  [coeff, FitInfo] = ridgeRegression(X, y, lambdaRatio, numCVFolds, numCVSamples, jitter) 

% [coeff, fitInfo] = ridgeRegression(X, y, lambdaRatop, numCVFolds, numCVSamples) 
% wrapper to perform ridge regression using a CV object and having similar
% output structure to lasso. Calls the matlab function ridge.m
% X is predictor matrix, y is response vector, lambdaRatio is vector
% containing the lambda hyperparameter. numCVfolds and numCVsamples are to
% set up k-fold CV object. alternatively numCVFolds can be the CV object
% itself
% coeff is n params x n lambda matrix for full data solution, fitInfo
% contains CV info similar to lasso and elasticNetRegression
%
% Original code found here: https://github.com/BrainCOGS/widefieldImaging

 %% Default arguments
if nargin < 3 || isempty(lambdaRatio)
  lambdaRatio           = flip(10.^(linspace(-3, log10(1), 30)'));
else
  lambdaRatio           = sort(lambdaRatio(:), 'descend');
end
if nargin < 4 || isempty(numCVFolds)
  numCVFolds            = 3;
end
if nargin < 5 || isempty(numCVSamples)
  numCVSamples          = 10;
end
if nargin < 6 || isempty(jitter)
  jitter                = [1e-6 1e-5];
end

%% Cross-validation setup
if isnumeric(numCVFolds)
  warning('off', 'stats:cvpartition:KFoldMissingGrp');
  pseudoExp           = cvpartition(y, 'KFold', numCVFolds);
  for iMC = 2:numCVSamples
    pseudoExp(iMC)    = repartition(pseudoExp(1));
  end
  warning('on', 'stats:cvpartition:KFoldMissingGrp');
else
  pseudoExp           = numCVFolds;
  numCVFolds          = pseudoExp(1).NumTestSets;
end

%% central value fits
FitInfo.CVExperiments     = pseudoExp;
FitInfo.Coefficients      = ridge(y, X, lambdaRatio);
coeff                     = FitInfo.Coefficients;

%% Cross validated goodness-of-fit
mse                       = nan(numel(lambdaRatio), numel(pseudoExp), numCVFolds);
rSquared                  = nan(numel(lambdaRatio), numel(pseudoExp), numCVFolds);
FitInfo.CVTrainPrediction = cell(numel(pseudoExp), numCVFolds);
FitInfo.CVTestPrediction  = cell(numel(pseudoExp), numCVFolds);

for iMC = 1:numel(pseudoExp)
  for iFold = 1:numCVFolds
    % Select train/test sets with user-specified pruning if so desired
    if isfield(pseudoExp(iMC).training(iFold),'idx')
      iTrain              = pseudoExp(iMC).training(iFold).idx;
      iTest               = pseudoExp(iMC).test(iFold).idx;
    else
      iTrain              = pseudoExp(iMC).training(iFold);
      iTest               = pseudoExp(iMC).test(iFold);
    end
    % Perform fit for this CV experiment and generate/test prediction
    mcCoeff               = ridge(y(iTrain), X(iTrain,:), lambdaRatio, 0); % zero also returns bias term

    % train prediction
    FitInfo.CVTrainPrediction{iMC,iFold} = ridgePrediction(X(iTrain,:),mcCoeff);
    [FitInfo.CVTestPrediction{iMC,iFold}, mse(:,iMC,iFold), rSquared(:,iMC,iFold)] ...
                          = ridgePrediction(X(iTest,:), mcCoeff, y(iTest), jitter);
  end
end

%% Collect output in same format as Matlab's lasso()
FitInfo.Lambda            = lambdaRatio;
FitInfo.DF                = sum(coeff ~= 0, 1);
FitInfo.MSE               = mean(mse(:,:), 2, 'omitnan')';
FitInfo.SE                = std(mse(:,:), 0, 2, 'omitnan')';
FitInfo.RSquared          = mean(rSquared(:,:), 2, 'omitnan')';
FitInfo.RSquaredErr       = std(rSquared(:,:), 0, 2, 'omitnan')';
[~, FitInfo.IndexMinMSE]  = min(FitInfo.MSE);
FitInfo.LambdaMinMSE      = FitInfo.Lambda(FitInfo.IndexMinMSE);
FitInfo.Index1SE          = find(   FitInfo.MSE(1:FitInfo.IndexMinMSE)    ...
                                <=  FitInfo.MSE(FitInfo.IndexMinMSE)      ...
                                  + FitInfo.SE(FitInfo.IndexMinMSE)       ...
                                , 1, 'first'                              ...
                                );
if isempty(FitInfo.Index1SE)
  FitInfo.Index1SE        = FitInfo.IndexMinMSE;
  FitInfo.Lambda1SE       = FitInfo.LambdaMinMSE;
else
  FitInfo.Lambda1SE       = FitInfo.Lambda(FitInfo.Index1SE);
end

end

%% predict data for different lambdas
function [yhat, mse, rSquared] = ridgePrediction(X, coeff, y, jitter)

if nargin < 3; y      = []; end
if nargin < 4; jitter = []; end

% add bias column
if size(X,2) < size(coeff,1)
  X = [ones(size(X,1),1) X];
end

nlambda = size(coeff,2);
npoints = size(X,1);
yhat    = zeros(npoints,nlambda);
for iL = 1:nlambda
  yhat(:,iL) = X*coeff(:,iL);
end

if isempty(y)
  mse      = [];
  rSquared = [];
else
  y        = repmat(y,[1 nlambda]);
  sse      = bsxfun(@minus, y, yhat).^2;
  sse      = sum(sse, 1);
  mse      = (sse ./ npoints)';
  muY      = repmat(mean(y),[npoints 1]);
  ssy      = sum((y - muY).^2);
  rSquared = (1 - sse ./ ssy)';
end

if ~isempty(jitter)
  yhat     = yhat + jitter(1)*randn(size(yhat));% + rand(size(yhat)) * (jitter(2) - jitter(1));
end

end