function [FREQ, selected_lambda] = randomized_elastic_net(expression,grp,nperm,nfolds,alpha,weakness,weakness2,lambda)
%% function [FREQ, selected_lambda] = randomized_lasso(expression,grp,nperm,nfolds,alpha,weakness,weakness2,lambda)
% expression genes x samples 
%               -> meaning observations in coulmns, genes in rows
% grp -> sample labels
% nperm -> number of permutations
% nfolds -> number of folds in cross-validation
% weakness -> for genes
% weakness -> for observations
% lambda-> optional


nsamples = length(grp);
xtmp = double(expression)';
% standard normalization of the expression:
xtmp  =  bsxfun(@times,bsxfun(@minus,xtmp,mean(xtmp)),1./std(xtmp));
% getting lambdas:
foldnumber = min([sum(grp==1),sum(grp==0)]);
if foldnumber<2
    foldnumber=2;
end
foldid =  crossvalind('Kfold', grp,foldnumber);
if isempty(nfolds)
	nfolds = min(length(grp),16);
end
nfolds= max(foldid);

%% %%% %%%%%%%%%%%%%%%%%%%%%%% %%% %%
% get lambda from cross-validation: %
% %%%% %%%%%%%%%%%%%%%%%%%%%%% %%%% %
if nargin<8
    [coef,CVerr] = lassoglm(xtmp,double(grp),'binomial','Alpha',alpha,'CV',nfolds);  
else
    [coef,CVerr] = lassoglm(xtmp,double(grp),'binomial','Lambda',lambda,'Alpha',alpha,'CV',nfolds);  
end
% lassoPlot(coef,CVerr,'plottype','CV');
selected_lambda = CVerr.Lambda1SE;

%% %%% %%%%%%%%%%%%%%%%%%%%%%%%%%% %%% %%
% randomized lasso for selected lambda: %
% %%%% %%%%%%%%%%%%%%%%%%%%%%%%%%% %%%% %

FREQ = zeros(size(xtmp,2),1);
parfor (p = 1 : nperm)
    W = rand(size(xtmp,2),1) .* (1-weakness)+ weakness;
    W= W/mean(W);
    R=  rand(size(xtmp,1),1) .* (1-weakness2)+ weakness2;
    R= R/mean(R);
    x= bsxfun(@times,xtmp,W');
    [coef,~] = lassoglm(x,double(grp),'binomial','Alpha',alpha,'Lambda',selected_lambda,'Weights',R,'Standardize',false);  
    FREQ = FREQ + double((coef~=0));   
end
