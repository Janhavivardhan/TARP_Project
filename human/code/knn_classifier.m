function prediction=knn_classifier(input,target,fe1)
%% KNN classifier
 %# target class 
trainClass = target;


%# compute pairwise distances between each test instance vs. all training data
D = pdist2(fe1', input, 'euclidean');
[D,idx] = sort(D, 2, 'ascend');

%# K nearest neighbors
K = 3;
D = D(:,1:K);
idx = idx(:,1:K);

%# majority vote
prediction = mode(trainClass(idx),2);
end
% disp(['Recognised Face  set =',num2str(prediction)]);