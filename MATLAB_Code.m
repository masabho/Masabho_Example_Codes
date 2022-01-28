%% Estimating HIV prevalence of Enumeration ID from Malawi MPHIA
% Author: Masabho Peter Milali
% Created: November 2021
% MATLAB R2019a

clear;
  clc;
   %fileName = 'Workspaces/PredictorVariablesFromMalawiHHSurvey_RemovedVIPScoreZeroAndLogicalPair.mat';
   fileName = 'Workspaces/Malawi_Ind_HH_CombinedVIPGreaterOrOne.mat';
   if (exist(fileName, 'file'))
       load(fileName);
   else 
      %[num,txt,raw] =
      %xlsread('PredictorVariablesFromMalawiHHSurvey_RemovedVIPScoreZeroAndLogicalPair.xlsx');%HHVariables
      [num,txt,raw] = xlsread('Malawi_Ind_HH_CombinedVIPGreaterOrOne.xlsx');
      save(fileName, 'num', 'txt', 'raw');
   end
   
%% Data preparation (Training and Testing sets)
%%% Training PLS binary classification model (threshold 0.1)

 DataForBinaryModel = num;
 C_ReplicatePrevalencePLS = 1; % setting number of replicates
 for iPrevalence = 1:C_ReplicatePrevalencePLS  
     SortedDataForBinaryModel = sortrows(DataForBinaryModel, 2); % sorted according to prevalence
     NumberColdSpots = length(SortedDataForBinaryModel(SortedDataForBinaryModel(:,2) < 0.05)) % Number of Enumeration Areas with Prevalences < 0.1
     EnuAreaColdSpots = SortedDataForBinaryModel(1: NumberColdSpots,:);
     EnuAreaColdSpots_New = EnuAreaColdSpots(:,3:end);
     EnuAreaColdSpots_New = [zeros(size((EnuAreaColdSpots_New),1),1),EnuAreaColdSpots_New];
     TrainSize_EnuAreaColdSpots_New = randsample(NumberColdSpots,(0.7*NumberColdSpots));
     Train_EnuAreaColdSpotsNew = EnuAreaColdSpots_New(ismember(1:NumberColdSpots, ... 
                                         TrainSize_EnuAreaColdSpots_New),:);% Training LessThanPointOne
     Test_EnuAreaColdSpotsNew = EnuAreaColdSpots_New(~ismember(1:NumberColdSpots, ...
                                         TrainSize_EnuAreaColdSpots_New),:); % Testing LessThanPointOne
     
     NumberHotSpots = length(SortedDataForBinaryModel(SortedDataForBinaryModel(:,2) >= 0.05));
     EnuAreaHotSpots = SortedDataForBinaryModel((NumberHotSpots+1): end,:);
     EnuAreaHotSpots_New = EnuAreaHotSpots(:,3:end);
     EnuAreaHotSpotsNew = [ones(size((EnuAreaHotSpots_New),1),1),EnuAreaHotSpots_New];
     TrainSize_EnuAreaHotSpotNew = randsample(NumberHotSpots,(0.7*NumberHotSpots));
     Train_EnuAreaHotSpotsNew = EnuAreaHotSpotsNew(ismember(1:NumberHotSpots, ... 
                                         TrainSize_EnuAreaHotSpotNew),:);% Training GreaterThanEqualPointOne
     Test_EnuAreaHotSpotsNew = EnuAreaHotSpotsNew(~ismember(1:NumberHotSpots, ...
                                         TrainSize_EnuAreaHotSpotNew),:); % Testing GreaterThanEqualPointOne
     
   %  PrevalenceCategorized = [EnuAreaPrevLessPointOne_New ; EnuAreaPrevGreaterEqualPointOne_New];
         
     Prevalence_train = [Train_EnuAreaColdSpotsNew; Train_EnuAreaHotSpotsNew];
     K_Class = randperm(size(Prevalence_train,1));
     RandomizedPrevalence_train = Prevalence_train(K_Class(1:size(Prevalence_train,1)),:);
     
     Prevalence_test = [Test_EnuAreaColdSpotsNew ; Test_EnuAreaHotSpotsNew];
     k_class = randperm(size(Prevalence_test,1));
     RandomizedPrevalence_test = Prevalence_test(k_class(1:size(Prevalence_test,1)),:);
     
     Prevalence_Xtrain_PLS = RandomizedPrevalence_train(:,2:end);
     Prevalence_Xtest_PLS = RandomizedPrevalence_test(:,2:end);
     Prevalence_Ytrain_PLS =  RandomizedPrevalence_train(:,1);
     Prevalence_Ytest_PLS = RandomizedPrevalence_test(:,1);
     
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%      %% Perform PLS regression with fifty components
%         [XtrainL, Ytrainl, XtrainS, YtrainS, MalawiPrevalence_beta, PCTVAR, MSE, stats] = ... 
%             plsregress(Prevalence_Xtrain_PLS, Prevalence_Ytrain_PLS, 50, 'CV', 10);
%     
%     %% Plot the percent of variance explained in the response variable as a function of the # of components
%     %http://www.mathworks.com/help/stats/examples/partial-least-squares-regression-and-principal-components-regression.html#zmw57dd0e5310
%         figure()
%         %subplot(2,1,1);
%         plot(1:50, cumsum(100*PCTVAR(2,:)), '-bo');
%         % RMSE = MSE(:,2:end);
%        % plot(1:30, cumsum(100*RMSE(2,:)), '-bo');
%         xlabel('Number of PLS components', 'FontSize', 16);
%         ylabel('% Variance Explained in the Output', 'FontSize', 16);
%         set(gca,'XTickLabel',[{'0'}, {'5'},{'10'},{'15'}, {'20'},{'25'},{'30'}],'FontSize', 14);
%         set(gcf, 'Position',[10,10,900, 1500]);
%         set(gca,'fontname','arial')
%           xlim([0 52]);
%           ylim([0 52]);
%           print(gcf,'PerformanceVsPLSComponents_MalawiInd&HHSurvey.tif','-dtiff','-r300');
%           
%         hold on % (plotting for An.arabiensis, calling from separate code)
%         
%         
%         %% Ploting regression coefficients weights against Variables
%          %Wavelength_Mayagaya = M_Mayagaya(:,1);
%          FeatureSpace = [3:47]';
%          col_size = size(FeatureSpace, size(MalawiPrevalence_beta,1)); 
%          NaN_row = NaN(1, col_size);
%         %FeatureSpace_Prev = [FeatureSpace; NaN_row];NaN row last
%          FeatureSpacePrev = [NaN_row; FeatureSpace];
%          figure()
%          subplot(2,1,1);
%          plot(FeatureSpacePrev,MalawiPrevalence_beta);
%          xlabel('FeatureSpace', 'FontSize',16);
%          ylabel('Regression Coefficients', 'FontSize', 16);
%          xlim([400 2400]);
%          set(gcf, 'Position',[10,10,900, 1500]);
%          set(gca,'fontname','arial')
%          
%          %% Calculate the normalized PLS weights.
%          % https://www.mathworks.com/help/stats/plsregress.html
%          W0 = stats.W ./ sqrt(sum(stats.W.^2,1));
%          % Normalizing beta coeffiecient like weights
%          NormalizedBeta = MalawiPrevalence_beta./sqrt(sum(MalawiPrevalence_beta.^2,1));
%          
%          % Calculate the VIP scores for ncomp components.
%          % https://www.mathworks.com/help/stats/plsregress.html
%          
%          p = size(XtrainL,1);
%          sumSq = sum(XtrainS.^2,1).*sum(Ytrainl.^2,1);
%          vipScore = sqrt(p* sum(sumSq.*(W0.^2),2) ./ sum(sumSq,2));
%          
%          % Finds variables with a VIP score greater than or eaqual to 1
%          % https://www.mathworks.com/help/stats/plsregress.html
%          indVIP = find(vipScore >= 1);
%          indvip = find(vipScore);
%          
%          % Plot the VIP scores.
%          figure()
%          scatter(1:length(vipScore),vipScore,'x')
%          hold on
%          scatter(indVIP,vipScore(indVIP),'rx')
%          plot([1 length(vipScore)],[1 1],'--k')
%          hold off
%          axis tight
%          xlabel('Position of Predictor Variables','FontSize', 16)
%          ylabel('VIP Scores', 'FontSize', 16)
%          set(gca,'XTickLabel',[{'0'}, {'20'},{'40'},{'60'}, {'80'},{'100'},{'120'},{'140'},{'160'}],'FontSize', 14);
%          set(gcf, 'Position',[10,10,900, 1500]);
%          print(gcf,'VIPScores_MalawiIndHHSurvey.tif','-dtiff','-r300');
%          
%          
     
 
     % Training classification model on ten component using 10-fold cross-validation
     [XtrainL, ytrainl, XtrainS, YtrainS, MalawiPrevalence_beta] = plsregress(Prevalence_Xtrain_PLS , Prevalence_Ytrain_PLS, 30, 'CV', 10);
   
     A_MalawiPrevalence = [ones(size(Prevalence_Xtest_PLS,1), 1),  Prevalence_Xtest_PLS];
     PredictedPrevalence_PLS = A_MalawiPrevalence * MalawiPrevalence_beta;
%      B_MalawiPrevalence = [ones(size(Prevalence_Xtrain_PLS,1),1), Prevalence_Xtrain_PLS];
%      PredictedPrevalence_Ytrain_PLS = B_MalawiPrevalence * MalawiPrevalence_beta;
    
     
     %Scoring the model accuracy
     
    Number_Ytest_ColdSpots_PLS(iPrevalence) = sum(Prevalence_Ytest_PLS == 0);
    Number_Ytest_HotSpots_PLS(iPrevalence) = sum(Prevalence_Ytest_PLS == 1);
     
    NumberCorrectlyPredicted_ColdSpots_PLS = sum((Prevalence_Ytest_PLS == 0)&(PredictedPrevalence_PLS < 0.5));
    NumberWronglyPredicted_ColdSpots_PLS = sum((Prevalence_Ytest_PLS == 1)&(PredictedPrevalence_PLS < 0.5));
    NumberCorrectlyPredicted_HotSpots_PLS = sum((Prevalence_Ytest_PLS == 1)&(PredictedPrevalence_PLS >= 0.5));
    NumberWronglyPredicted_HotSpots_PLS = sum((Prevalence_Ytest_PLS == 0)&(PredictedPrevalence_PLS >= 0.5));
    AccuracyPrevalenceEstimate_PLS = ((NumberCorrectlyPredicted_ColdSpots_PLS + NumberCorrectlyPredicted_HotSpots_PLS)/length(Prevalence_Ytest_PLS));
    
    TrueColdSpots_PLS(iPrevalence) = NumberCorrectlyPredicted_ColdSpots_PLS;
    FalseColdSpots_PLS(iPrevalence) = NumberWronglyPredicted_ColdSpots_PLS;
    TrueHotSpots_PLS(iPrevalence) =  NumberCorrectlyPredicted_HotSpots_PLS;
    FalseHotSpots_PLS(iPrevalence) = NumberWronglyPredicted_HotSpots_PLS;
    AccuracyPrevalenceEstimatePLS(iPrevalence) = AccuracyPrevalenceEstimate_PLS;     
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %% Average Scores of PLS model on estimating Prevalence     
     AverageNumber_Ytest_ColdSpots_PLS = mean(Number_Ytest_ColdSpots_PLS);
     NumberPredicted_ColdSpots_PLS =  TrueColdSpots_PLS  + FalseColdSpots_PLS;
     AverageNumberPredicted_ColdSpots_PLS = mean(NumberPredicted_ColdSpots_PLS);
     AverageNumber_Ytest_HotSpots_PLS = mean(Number_Ytest_HotSpots_PLS);
     NumberPredicted_HotSpots_PLS = TrueHotSpots_PLS + FalseHotSpots_PLS;
     AverageNumberPredicted_HotSpots_PLS = mean(NumberPredicted_HotSpots_PLS);
     
     Average_Accuracy_Prevalence_PLS = mean(AccuracyPrevalenceEstimatePLS);
     Std_dev_Accuracy_Prevalence_PLS = std(AccuracyPrevalenceEstimatePLS)
     
     Average_TrueColdSpots_PLS = mean(TrueColdSpots_PLS);
     Average_FalseColdSpots_PLS = mean(FalseColdSpots_PLS);
     Average_TrueHotSpots_PLS =  mean(TrueHotSpots_PLS);
     Average_FalseHotSpots_PLS = mean(FalseHotSpots_PLS);
      
     Sensitivity_Prevalence_PLS = TrueHotSpots_PLS./Number_Ytest_HotSpots_PLS
     Average_Sensitivity_Prevalence_PLS = mean(Sensitivity_Prevalence_PLS)
     StdDev_Sensitivity_Prevalence_PLS = std(Sensitivity_Prevalence_PLS)
     
     Specificity_Prevalence_PLS = TrueColdSpots_PLS./Number_Ytest_ColdSpots_PLS
     Average_Specificity_Prevalence_PLS = mean(Specificity_Prevalence_PLS)
     StdDev_Specificity_Prevalence_PLS = std(Specificity_Prevalence_PLS)
     
     Precision_Prevalence_PLS = TrueHotSpots_PLS./(TrueHotSpots_PLS + FalseHotSpots_PLS)
     Average_Precision_Prevalence_PLS = mean(Precision_Prevalence_PLS)
     stdDev_Precision_Prevalence_PLS = std(Precision_Prevalence_PLS)
     
     Recall_Prevalence_PLS = TrueHotSpots_PLS./(TrueHotSpots_PLS + FalseColdSpots_PLS)
     Average_Recall_Prevalence_PLS = mean(Recall_Prevalence_PLS)
     stdDev_Recall_Prevalence_PLS = std(Recall_Prevalence_PLS)
%      
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   %% Plotting Figures 
%     
%     %% box plot 
%     figure()
%     set(gcf, 'Position',[10,10,850, 2000]);
%     subplot(2,2,1);
%     %subplot(3,2,4);
%     h = boxplot(PredictedPrevalence_PLS, Prevalence_Ytest_PLS, 'Colors','k');
%     set(h,{'linew'},{2})
%     ylim([-0.3 1.3])
%     xlabel('Actual Prevalence', 'FontSize',14);
%     ylabel('Estimated Prevalence', 'FontSize', 14);
%     set(gca,'XTickLabel',[{'ColdSpot (0)'}, {'HotSpot (1)'}],'FontSize', 14);
% %     set(gca,'YTickLabel',[{''}, {'-0.5'}, {'0'},{'0.5'},{'1'},{''}],'FontSize', 14);
%     set(gca,'fontname','arial');
% %     print(gcf,'BoxPlotANN_HIVPrevalenceFig_.tif','-dtiff','-r300');
% 
% 
%      %%  bar graph
%      figure()
%      %subplot(2,2,2);
%     % subplot(2,2,3);
%      xdata_Prevalence = [1 2];
%      ydata_Prevalence = [Average_TrueColdSpots_PLS/AverageNumberPredicted_ColdSpots_PLS ... 
%          *100 Average_FalseColdSpots_PLS/AverageNumberPredicted_ColdSpots_PLS*100; ...
%          Average_TrueHotSpots_PLS/AverageNumberPredicted_HotSpots_PLS*100 ... 
%          Average_FalseHotSpots_PLS/AverageNumberPredicted_HotSpots_PLS*100];
%      E = [std(TrueColdSpots_PLS./NumberPredicted_ColdSpots_PLS)*100, ... 
%          std(FalseColdSpots_PLS./NumberPredicted_ColdSpots_PLS)*100; ... 
%          std(TrueHotSpots_PLS./NumberPredicted_HotSpots_PLS)*100, ... 
%          std(FalseHotSpots_PLS./NumberPredicted_HotSpots_PLS)*100];
%      hb = bar(xdata_Prevalence,ydata_Prevalence,1);
%      hold on
%      xe = [0.85 1.15; 1.85 2.15];
%      errorbar(xe,ydata_Prevalence,E,'*','CapSize',17,'LineWidth',2)
%      
%      set(hb(2),'facecolor',[244/256 134/256 66/256])
%      %set(hb(1),'facecolor',[179/256 66/256  244/256])
%      %set(hb(1),'facecolor',[65/256 244/256 244/256])
%      set(hb(1),'facecolor',[176/256 244/256 66/256])
%      ylabel('% Estimated prevalence','FontSize', 14);
%      xlabel('Model prediction','FontSize', 14); 
%      %set(gca,'XTick',[1 2]);
%      set(gca,'XTickLabel',[{'ColdSpots (N = 76)'}, {'HotSpots (N = 75)'}],'FontSize', 14);
%      ylim([0 140])
%      set(gca,'YTickLabel',[{'0'}, {'20'},{'40'},{'60'}, {'80'},{'100'},{''},{''}],'FontSize', 14);
%      
%      legend({'Correct prediction','False prediction'},'FontSize',14, ...
%               'location','NW');
%      legend boxoff
%      set(gca,'fontname','arial')
%      print(gcf,'BarPlot_Prevalence_MalawiHHSurvey_PLS_.tif','-dtiff','-r300');
%      hold off
%      
% %% ROC curve
% % https://stackoverflow.com/questions/33523931/matlab-generate-confusion-matrix-from-classifier/33542453#33542453
%      tot_op =  PredictedPrevalence_PLS;
%      targets = Prevalence_Ytest_PLS;
%      th_vals= sort(tot_op);
% 
%      for i = 1:length(th_vals)
%        b_pred = (tot_op>=th_vals(i,1));
%        TP = sum(b_pred >=0.5 & targets == 1);
%        FP = sum(b_pred >=0.5 & targets == 0);
%        TN = sum(b_pred < 0.5 & targets == 0);
%        FN = sum(b_pred < 0.5 & targets == 1);
%         sens(i) = TP/(TP+FN);
%         spec(i) = TN/(TN+FP);
%      end
% 
% 
%       figure()
%       %subplot(2,2,4);
%       cspec = 1-spec;
%       cspec = cspec(end:-1:1);
%       sens = sens(end:-1:1);
%       plot(cspec,sens, 'k','linew',2);
%       hold on
%       %legend('Minepa-ARA','Muleba-GA')
%       xlabel('False positive rate','FontSize',14); ylabel('True positive rate','FontSize',14);
% %       set(gca,'YTickLabel',[{'0'}, {'0.2'}, {'0.4'},{'0.6'},{'0.8'},{'1'}],'FontSize', 14);
% %       title('ROC Curves for Minepa-GA, Muleba-GA, and Burkina-GA and Muleba-Burkina-GA')
%      % print(gcf,'FieldPaperROC_CURVE.tif','-dtiff','-r300');
%       % print(gcf,'FieldPaperROC_CURVE_Autoencoder_Purtubed.tif','-dtiff','-r300');
%       print(gcf,'AOC_Prevalence_PLS_.tif','-dtiff','-r300');
% 
%       AUC = sum(0.5*(sens(2:end)+sens(1:end-1)).*(cspec(2:end) - cspec(1:end-1)));
%       fprintf('\nAUC: %g \n',AUC);
      
      
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Training ANN
  DataForBinaryModel = num;
 C_ReplicatePrevalenceANN = 10; % setting number of replicates
 for iPrevalence = 1:C_ReplicatePrevalenceANN  
     SortedDataForBinaryModel = sortrows(DataForBinaryModel, 2); % sorted according to prevalence
     NumberColdSpots = length(SortedDataForBinaryModel(SortedDataForBinaryModel(:,2) < 0.1)) % Number of Enumeration Areas with Prevalences < 0.1
     EnuAreaColdSpots = SortedDataForBinaryModel(1: NumberColdSpots,:);
     EnuAreaColdSpots_New = EnuAreaColdSpots(:,3:end);
     EnuAreaColdSpots_New = [zeros(size((EnuAreaColdSpots_New),1),1),EnuAreaColdSpots_New];
     TrainSize_EnuAreaColdSpots_New = randsample(NumberColdSpots,(0.7*NumberColdSpots));
     Train_EnuAreaColdSpotsNew = EnuAreaColdSpots_New(ismember(1:NumberColdSpots, ... 
                                         TrainSize_EnuAreaColdSpots_New),:);% Training LessThanPointOne
     Test_EnuAreaColdSpotsNew = EnuAreaColdSpots_New(~ismember(1:NumberColdSpots, ...
                                         TrainSize_EnuAreaColdSpots_New),:); % Testing LessThanPointOne
     
     NumberHotSpots = length(SortedDataForBinaryModel(SortedDataForBinaryModel(:,2) >= 0.1));
     EnuAreaHotSpots = SortedDataForBinaryModel((NumberHotSpots+1): end,:);
     EnuAreaHotSpots_New = EnuAreaHotSpots(:,3:end);
     EnuAreaHotSpotsNew = [ones(size((EnuAreaHotSpots_New),1),1),EnuAreaHotSpots_New];
     TrainSize_EnuAreaHotSpotNew = randsample(NumberHotSpots,(0.7*NumberHotSpots));
     Train_EnuAreaHotSpotsNew = EnuAreaHotSpotsNew(ismember(1:NumberHotSpots, ... 
                                         TrainSize_EnuAreaHotSpotNew),:);% Training GreaterThanEqualPointOne
     Test_EnuAreaHotSpotsNew = EnuAreaHotSpotsNew(~ismember(1:NumberHotSpots, ...
                                         TrainSize_EnuAreaHotSpotNew),:); % Testing GreaterThanEqualPointOne
     
   %  PrevalenceCategorized = [EnuAreaPrevLessPointOne_New ; EnuAreaPrevGreaterEqualPointOne_New];
         
     Prevalence_train = [Train_EnuAreaColdSpotsNew; Train_EnuAreaHotSpotsNew];
     K_Class = randperm(size(Prevalence_train,1));
     RandomizedPrevalence_train = Prevalence_train(K_Class(1:size(Prevalence_train,1)),:);
     
     Prevalence_test = [Test_EnuAreaColdSpotsNew ; Test_EnuAreaHotSpotsNew];
     k_class = randperm(size(Prevalence_test,1));
     RandomizedPrevalence_test = Prevalence_test(k_class(1:size(Prevalence_test,1)),:);
     
     Prevalence_Xtrain_ANN = RandomizedPrevalence_train(:,2:end);
     Prevalence_Xtest_ANN = RandomizedPrevalence_test(:,2:end);
     Prevalence_Ytrain_ANN =  RandomizedPrevalence_train(:,1);
     Prevalence_Ytest_ANN = RandomizedPrevalence_test(:,1);
     
     netPrevalenceClass = feedforwardnet(10);       % One hidden layer with ten neurons
%     netParityStatus.layers{end}.transferFcn = 'logsig';    % logistic regression as transfer function
    %netParityStatus.performFcn = 'crossentropy';        % default matlab uses Levenberg-Marquardt optimization
    %[trainInd,valInd,testInd] = dividerand(637,0.7,0.15,0.15); % customizing size of training, validation and test sets
    

   %% Training of the neural network model
   [netPrevalenceClass,tr] = train(netPrevalenceClass,Prevalence_Xtrain_ANN',Prevalence_Ytrain_ANN');
   
    Prevalence_TrainPredicted = netPrevalenceClass(Prevalence_Xtrain_ANN');% training dataset
    PredictedPrevalence_ANN = netPrevalenceClass(Prevalence_Xtest_ANN');%out of the sample data set
    PredictedPrevalence_ANN = PredictedPrevalence_ANN';

   %% Scoring the accuracy of ANN Prevalence Classification Model
    
    Number_Ytest_ColdSpots_ANN(iPrevalence) = sum(Prevalence_Ytest_ANN == 0);
    Number_Ytest_HotSpots_ANN(iPrevalence) = sum(Prevalence_Ytest_ANN == 1);
     
    NumberCorrectlyPredicted_ColdSpots_ANN = sum((Prevalence_Ytest_ANN == 0)&(PredictedPrevalence_ANN < 0.5));
    NumberWronglyPredicted_ColdSpots_ANN = sum((Prevalence_Ytest_ANN == 1)&(PredictedPrevalence_ANN < 0.5));
    NumberCorrectlyPredicted_HotSpots_ANN = sum((Prevalence_Ytest_ANN == 1)&(PredictedPrevalence_ANN >= 0.5));
    NumberWronglyPredicted_HotSpots_ANN = sum((Prevalence_Ytest_ANN == 0)&(PredictedPrevalence_ANN >= 0.5));
    AccuracyPrevalenceEstimate_ANN = ((NumberCorrectlyPredicted_ColdSpots_ANN + NumberCorrectlyPredicted_HotSpots_ANN)/length(Prevalence_Ytest_ANN));
    
    TrueColdSpots_ANN(iPrevalence) = NumberCorrectlyPredicted_ColdSpots_ANN;
    FalseColdSpots_ANN(iPrevalence) = NumberWronglyPredicted_ColdSpots_ANN;
    TrueHotSpots_ANN(iPrevalence) =  NumberCorrectlyPredicted_HotSpots_ANN;
    FalseHotSpots_ANN(iPrevalence) = NumberWronglyPredicted_HotSpots_ANN;
    AccuracyPrevalenceEstimateANN(iPrevalence) = AccuracyPrevalenceEstimate_ANN;
    
 end 
  
 
 %% Average Scores of PLS model on estimating Prevalence     
     AverageNumber_Ytest_ColdSpots_ANN = mean(Number_Ytest_ColdSpots_ANN);
     NumberPredicted_ColdSpots_ANN =  TrueColdSpots_ANN  + FalseColdSpots_ANN;
     AverageNumberPredicted_ColdSpots_ANN = mean(NumberPredicted_ColdSpots_ANN);
     AverageNumber_Ytest_HotSpots_ANN = mean(Number_Ytest_HotSpots_ANN);
     NumberPredicted_HotSpots_ANN = TrueHotSpots_PLS + FalseHotSpots_ANN;
     AverageNumberPredicted_HotSpots_ANN = mean(NumberPredicted_HotSpots_ANN);
     
     Average_Accuracy_Prevalence_ANN = mean(AccuracyPrevalenceEstimateANN);
     Std_dev_Accuracy_Prevalence_ANN = std(AccuracyPrevalenceEstimateANN)
     
     Average_TrueColdSpots_ANN = mean(TrueColdSpots_ANN);
     Average_FalseColdSpots_ANN = mean(FalseColdSpots_ANN);
     Average_TrueHotSpots_ANN =  mean(TrueHotSpots_ANN);
     Average_FalseHotSpots_ANN = mean(FalseHotSpots_ANN);
      
     Sensitivity_Prevalence_ANN = TrueHotSpots_PLS./Number_Ytest_HotSpots_ANN
     Average_Sensitivity_Prevalence_ANN = mean(Sensitivity_Prevalence_ANN)
     StdDev_Sensitivity_Prevalence_ANN = std(Sensitivity_Prevalence_ANN)
     
     Specificity_Prevalence_ANN = TrueColdSpots_PLS./Number_Ytest_ColdSpots_ANN
     Average_Specificity_Prevalence_ANN = mean(Specificity_Prevalence_ANN)
     StdDev_Specificity_Prevalence_ANN = std(Specificity_Prevalence_ANN)
     
     Precision_Prevalence_ANN = TrueHotSpots_ANN./(TrueHotSpots_ANN + FalseHotSpots_ANN)
     Average_Precision_Prevalence_ANN = mean(Precision_Prevalence_ANN)
     stdDev_Precision_Prevalence_ANN = std(Precision_Prevalence_ANN)
     
     Recall_Prevalence_ANN = TrueHotSpots_ANN./(TrueHotSpots_ANN + FalseColdSpots_ANN)
     Average_Recall_Prevalence_ANN = mean(Recall_Prevalence_ANN)
     stdDev_Recall_Prevalence_ANN = std(Recall_Prevalence_ANN)
