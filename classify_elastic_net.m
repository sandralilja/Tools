function [AUC,PVAL] = classify_elastic_net(xtmp, grp,selected_lambda)
% xtmp - expression matrix genes in rows, samples in columns
% grp - sample labels
% selected_lambda - selected lambda value

foldnumber = min([sum(grp==1),sum(grp==0)]);
foldid =  crossvalind('Kfold', grp,foldnumber);
nsamples= length(grp);
nfolds= max(foldid);
label  = zeros(nsamples,1);
xtmp = xtmp';
xtmp  =  bsxfun(@times,bsxfun(@minus,xtmp,mean(xtmp)),1./std(xtmp));
for k=1:nfolds;
    incv = true(nsamples,1);
    incv(foldid==k) = false;
    [coefN,CVerrN] = lassoglm(xtmp(incv,:),double(grp(incv)),'binomial','Alpha',0.5,'CV',nfolds,'Lambda',selected_lambda);
    tmp  =find(~incv);
    for j=tmp
        label(j) = glmval([CVerrN.Intercept(CVerrN.IndexMinDeviance);coefN(:,CVerrN.IndexMinDeviance)],xtmp(~incv,:),'logit');
    end
end
[~,~,~,AUC]= perfcurve(logical(grp),label,true);
PVAL = ranksum( label(grp ==1),label(grp==0));

