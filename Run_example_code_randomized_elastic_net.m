clear all
close all
clc

%% read in example data set
% This is breast cancer data we have used for validation of the classifier SHARED project
% genes were filetered by the presence in the BC disease module

data.table = readtable('exampleData.txt');
data.data= table2array(data.table(:,2:end)); % becouse first column contains entrez gene IDs
data.class = table(data.table.Properties.VariableNames(2:end)','variablenames',{'title'});
data.class.BC = ~cellfun(@isempty,strfind(data.class.title,'BC'));% logical, frue for BC patients, false for controls, identifed by search of the 'BC' string in column names


%% run randomized elastic net (if you have many genes this might take a while...)
% gene prefiltering might be needed like doing it only for module genes,
% or for genes predicted to be secreted (Human Protein Atlas) or something
% else that is relevant in your project ;)
% The example expression data contains only breast cancer module genes.

help randomized_lasso % function description

% now I will exaggerate:
expression = data.data; % data matrix samples in columns, genes in rows
grp = data.class.BC; % classes
nperm = 1000; % number of permutations that is what we used in SHARED ????
nfolds =16; % how many folds do you want to have? max is 16 it's in the function code so you can change the max whenever you want ;)
alpha =0.5; % shift between ridge and lasso - this value we have used in SHARED and other analyses
weakness = 0.1; % that we have used in SHARED
weakness2 = 0.1; % that we have used in SHARED
% if you want you can also specify lambda instead of choosing it from
% cross-validation - I skip it.
[FREQ, selected_lambda] = randomized_elastic_net(expression,grp,nperm,nfolds,alpha,weakness,weakness2)
Results = table(data.table.probe, FREQ./nperm , 'variablenames',{'GeneName','frequency'});
writetable(Results,'ExampleResults.txt','delimiter','\t')

%% if you want to continue to do classifications:
% here you won't get exactly same results as reported in SHARED manuscript 
% because there this data was classified using lambda and genes selected
% based on different BC patients and healthy control cohort.
expression(FREQ==0,:)=[];
[AUC,PVAL] = classify_elastic_net(expression, grp,selected_lambda)

