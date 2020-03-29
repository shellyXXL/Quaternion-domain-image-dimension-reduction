function [ cor_ratio ] = slda_img_2d(persons,number, train_number,name,trainindex, testindex, kRange,stopNum, miu)
%PCA_IMG Summary of this function goes here
%   Detailed explanation goes here

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

class_mean=[];
for rr=1:persons
    start_ind=(rr-1)*train_number+1;
    end_ind=rr*train_number;
    tempt_bm=mean(sample_quater(start_ind:end_ind,:),1);
    class_mean=[class_mean; tempt_bm];
end
sample_mean=mean(sample_quater,1); 
%%%%%%%%%%%% calculate sw  sb  from  left
t=zeros(size_a,size_a);
sb=quaternion(t,t,t,t);
sw=quaternion(t,t,t,t);
for rr=1:persons
    tempS=quaternion(t,t,t,t);
    for ss=1:train_number
         img_1d=sample_quater((rr-1)*train_number+ss,:)-class_mean(rr,:);
         img_2d=reshape(img_1d, size_a, size_b);
         tempS=tempS+img_2d*img_2d';    
    end
    sw=sw+tempS./train_number;
    img_1d=class_mean(rr,:)-sample_mean;
    img_2d=reshape(img_1d, size_a, size_b);
    sb=sb+train_number*img_2d*img_2d'; 
end
sw=sw./persons;
sb=sb./(persons*train_number);

p=eye(size_a,size_a);
pp=quaternion(p,t,t,t);
sw=sw+0.0001*pp;


%% main procedure 
delta = 0;
stops_AB= 0.001;
maxSteps_AB = 5;
maxSteps_admm=100;

[ V] = qslda(sb, sw, K, delta, maxSteps_AB, stops_AB, maxSteps_admm,stopNum, miu);


  %% testing data 
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

clear class_mean;
clear sb;
clear sw;
clear img_2d;
clear cov;



%%%%%%%%%% classification
for y=1:K/kRange
    warning('off','all') 
    k=kRange*y;    
     bases=V(:,1:k);

    %% extract training features
    t=zeros(persons,train_number,k,size_b);
    train_feature=quaternion(t,t,t,t);
    for rr=1:persons
        for ss=1:train_number
             ind=(rr-1)*train_number+ss;
            temp_train_2d=reshape(sample_quater(ind,:), size_a, size_b);
                 result=bases'*temp_train_2d;
             train_feature(rr,ss,:,:)=result;     
        end
    end

    %% extract testing features 
    t=zeros(persons,test_number,k,size_b);
    test_feature=quaternion(t,t,t,t);
      for rr=1:persons
        for ss=1:test_number
             ind=(rr-1)*test_number+ss;
            temp_test_2d=reshape(test_quater(ind,:), size_a, size_b);
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
            
                 [value0, index0]=min(distance0);
                 distance1(mm)=value0;
            end
                     
            [value1, index1]=min(distance1);
            if index1~=rr
                errors=errors+1;

            end
        end
    end
    
    cor_ratio(y)=(persons*test_number-errors)/(persons*test_number)   
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

