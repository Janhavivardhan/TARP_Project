function [b] = svm_classifier(TrainingSet,GroupTrain,TestSet,sel)
%Models a given training set with a corresponding group vector and 
%classifies a given test set using an SVM classifier 

u=unique(GroupTrain);
numClasses=length(u);
result = zeros(length(TestSet(1,:)),1);
%%
% disp('Training The Data ....');
model = cell(length(u),1);
for k=1:numClasses
    a=double(GroupTrain==k);
    
%     a(find(a==0))=2;
    model{k} = svmtrain( TrainingSet,a,'kernel_function','rbf','RBF_Sigma', 0.6,'BoxConstraint', 0.8);
end
out=zeros(1,length(GroupTrain))';
%# get probability estimates of test instances using each model
% prob = zeros(numTest,numLabels); 
ss=0;
% disp('Testing The Data ....');
if sel==1
for k=1:numClasses
    p(k) = svmclassify(model{k},TestSet');
    
%     prob(:,k) = p(:,model{k}.Label==1);    %# probability of class==k
end
%% 
[a,b]=max(p);
else
   for jk=1:size(TestSet,2) 
    for k=1:numClasses
    p(k) = svmclassify(model{k},TestSet(:,jk)');
    
%     prob(:,k) = p(:,model{k}.Label==1);    %# probability of class==k
    end
%% 

[a,b1]=max(p);
b(jk)=b1;
end
end

