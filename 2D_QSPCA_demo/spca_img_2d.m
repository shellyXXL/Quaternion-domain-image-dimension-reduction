function [ cor_ratio ] = spca_img_2d(persons,number, train_number,name,stopNum, trainindex, testindex, kRange)

trainInd1=cell2mat(trainindex(1));
trainInd2=cell2mat(trainindex(2));
testInd1=cell2mat(testindex(1));
testInd2=cell2mat(testindex(2));

%%%%%% normalization data to 0-1
load([name '.mat']);
if mean(mean(all_sample.x))>1
   all_sample=all_sample./255;
end
size_a=data_size(1);
size_b=data_size(2);
K=size_a; %%%%%%%%%% number of projection vectors

sample_quater=[];
for gg=1:persons
    for hhh=trainInd1(1):trainInd1(2)
         ind=(gg-1)*number+hhh;
         temp=all_sample(ind,:);
         sample_quater=[sample_quater ;temp];
    end
    if ~isempty(trainInd2)
       for hhh=trainInd2(1):trainInd2(2)
         ind=(gg-1)*number+hhh;
         temp=all_sample(ind,:);
         sample_quater=[sample_quater ; temp];
       end
    end
end

%%normalization
sample_mean=mean(sample_quater,1);
trainTatal=persons*train_number;
sample_meanM=repmat(sample_mean,trainTatal,1);
sample_quaterAdj=sample_quater-sample_meanM;

%%2d image processing
nn=persons*train_number;
tmp_size=zeros(size_a);
cov=quaternion(tmp_size,tmp_size,tmp_size,tmp_size);
for i=1:nn
    img_2d=reshape(sample_quaterAdj(i,:), size_a, size_b);
    cov=cov+img_2d*img_2d';
end
cov=cov./nn;

%%main sparse QSPCA procedure
delta = 0;
stops_AB= 0.001;
maxSteps_AB = 20;
maxSteps_admm=100;   
ro=0.001;
Vs = qspca(cov, K, delta, maxSteps_AB, stops_AB, ro,maxSteps_admm,stopNum);

%%extract testing data 
test_number=testInd1(2)-testInd1(1)+1;
 if ~isempty(testInd2)
     test_number=test_number+testInd2(2)-testInd2(1)+1;
 end
test_quater=[];
for gg=1:persons
 for hhh=testInd1(1):testInd1(2)
         ind=(gg-1)*number+hhh;
         temp=all_sample(ind,:);
        test_quater=[test_quater; temp];
 end   
    if ~isempty(testInd2)
       for hhh=testInd2(1):testInd2(2)
         ind=(gg-1)*number+hhh;
         temp=all_sample(ind,:);
        test_quater=[test_quater;temp];
       end
    end
end

%%normalization testing data
testTatal=persons*test_number;
test_meanM=repmat(sample_mean,testTatal,1);
test_quaterAdj=test_quater-test_meanM;

clear all_sample;
clear test_meanM;
clear sample_meanM;
clear sample_quater;
clear test_quater;
clear img_2d;
clear cov;
clear tmp_size;

for y=1:K/kRange 
    k=kRange*y;
    bases=Vs(:,1:k);
    %%extract training features
    t=zeros(persons,train_number,k,size_b);
    train_feature=quaternion(t,t,t,t);
    for rr=1:persons
        for ss=1:train_number
             ind=(rr-1)*train_number+ss;
            temp_train_2d=reshape(sample_quaterAdj(ind,:), size_a, size_b);
                 result=bases'*temp_train_2d;
             train_feature(rr,ss,:,:)=result;     
        end
    end

    %%extract testing features 
    t=zeros(persons,test_number,k,size_b);
    test_feature=quaternion(t,t,t,t);
    for rr=1:persons
        for ss=1:test_number
             ind=(rr-1)*test_number+ss;
            temp_test_2d=reshape(test_quaterAdj(ind,:), size_a, size_b);
             result=bases'*temp_test_2d;
             test_feature(rr,ss,:,:)=result;
        end
    end
     %%% nearest neighbor classifer.
    errors=0;
    for rr=1:persons
        for ss=1:test_number
            my_tempt2=test_feature(rr,ss,:,:);
            distance1=zeros(1,persons);
            for mm=1:persons
                distance0=zeros(1,train_number); 
                for nn=1:train_number
                    my_tempt1=train_feature(mm,nn,:,:);
                    dis_q=(my_tempt2(:)-my_tempt1(:));
                    distance0(nn)=norm(dis_q,1);
                end
                [value0, ~]=min(distance0);
                 distance1(mm)=value0;
            end
            
            [~, index1]=min(distance1);
            if index1~=rr
                errors=errors+1;
            end           
        end      
    end
    
    cor_ratio(y)=(persons*test_number-errors)/(persons*test_number);   
        
end
    clear train_feature;
    clear test_feature;
    clear my_tempt2;
    clear my_tempt1;
    clear  distance1;
    clear  distance0;
    clear dis_q;
    clear temp_test_2d;
    clear temp_train_2d;
    clear result;

end

