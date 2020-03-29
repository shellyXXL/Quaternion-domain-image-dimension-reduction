%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  demo of "Two-dimensional quaternion sparse discriminant analysis", TIP 2018.
%%%%%  by Xiaolin Xiao, email:shellyxiaolin@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc, 
close all
warning off;
addpath(genpath('qtfm'));

kRange=2; %%%%%%%%%%%%%%%%%%%%% test classification accuracy per kRange for efficiency, you can also set kRange to 1.
stopNumArr=2:2:16;   % number of nonzero loadings per basis vector, usually, very sparse vectors are enough, if not you can test all possible solutions, e.g., stopNumArr=2:2:32
miuArr=[10e-3, 10e-2,10e-1];   % ratio of Sb Sw, usually [10e-3, 10e-2,10e-1] is enough, if not, you can test  [10e-3, 10e-2,10e-1,10e0, 10e1,10e2,10e3]
%%%%%%%%%%%%%%%%%%%%%%%%%%  stopNumArr=2;  
number=14;
persons=52;    
dataset='EU_cd';  %%%%%%%%%%%%  RGBD images where depth in real part  

%%%% train and test index, on clean EU 
trainindex={[1,4],[]};  testindex={[8,11],[]};   %%%%%% best parameter stopNum=6, miu=10e-2
% % % %%%% train and test index, on occluded EU  
% % % trainindex={[1,4],[8,11]};  testindex={[5,7],[12,14]};   %%%%%% best parameter stopNum=2, miu=10e-2
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
for j=1:length(miuArr)
     miu= miuArr(j); 
    for i=1:length(stopNumArr)  
       stopNum=stopNumArr(i); 
        cor_ratio_v = slda_img_2d(persons,number, train_number,dataset,trainindex, testindex, kRange,stopNum, miu);
        result_v=[result_v;cor_ratio_v];
    end
end
saveName=['2DQSDA_'   dataset   '_train' num2str(train_number) '_test' num2str(test_number) '.txt'];
dlmwrite(saveName, result_v);




