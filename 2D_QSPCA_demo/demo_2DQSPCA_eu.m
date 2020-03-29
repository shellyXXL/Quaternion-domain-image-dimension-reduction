%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  demo of "Two-dimensional quaternion PCA and sparse PCA", TNNLS 2018.
%%%%%  by Xiaolin Xiao, email:shellyxiaolin@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc, close all
warning off;
addpath(genpath('qtfm'));

kRange=2; %%%%%%%%%%%%%%%%%%%%% test classification accuracy per kRange for efficiency, you can also set kRange to 1.
stopNumArr=2:2:20;   %%%%%%%%%%%%%%%%%%%%% number of nonzero loadings per basis vector  full test stopNumArr=2:2:32;  
dataset='EU';   %%%%%%%%%% 14 images per subject, 52 subjects
number=14;
persons=52;
%%%% train and test index, on clean EU 
trainindex={[1,4],[]};  testindex={[8,11],[]}; 
% % %%%% train and test index, on occluded EU 
% % trainindex={[1,4],[8,11]};  testindex={[5,7],[12,14]}; 
trainInd1=cell2mat(trainindex(1));
trainInd2=cell2mat(trainindex(2));
testInd1=cell2mat(testindex(1));
testInd2=cell2mat(testindex(2));
train_number=trainInd1(2)-trainInd1(1)+1;
 if ~isempty(trainInd2)
     train_number=train_number+trainInd2(2)-trainInd2(1)+1;
 end
  test_number=testInd1(2)-testInd1(1)+1;
if ~isempty(testInd2)
     test_number=test_number+testInd2(2)-testInd2(1)+1;
 end

result_v=[];
for i=1:length(stopNumArr)
     stopNum=stopNumArr(i);
    [cor_ratio_v] = spca_img_2d(persons,number, train_number,dataset,stopNum, trainindex, testindex, kRange)
    result_v=[result_v;cor_ratio_v];
end
saveName=['2DQSPCA_'   dataset   '_train' num2str(train_number) '_test' num2str(test_number) '.txt'];
dlmwrite(saveName, result_v);




