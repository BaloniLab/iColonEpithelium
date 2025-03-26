%% GEM per sample UC
% healthy sample
folderPath = 'C:\Users\Alvis\Box\Boyu-Gut_microbiome\CyberGut\Case_Studies\Broad_Study\FAV_analysis_per_sample\Rxn_scores_0518_2024\Healthy_RxnScore\'
files = dir(folderPath)
fileNames = {files(~[files.isdir]).name};
model = readCbModel('C:\Users\Alvis\Box\Boyu-Gut_microbiome\CyberGut\Draft_reconstructions\Draft_reconstruction_consensus\Upadate_05-08-2024\icolonoEpithelium_loopless.mat');
DMEM = table2cell(DMEMinputs);
model_EXCHANGE = strmatch('EX_', model.rxns);
model.lb(find(ismember(model.rxns,model.rxns(model_EXCHANGE))))=0;
model.ub(find(ismember(model.rxns,model.rxns(model_EXCHANGE))))=1000;

for i = 1:1:58
    dmemreaction = DMEM{i}
    model.lb(find(ismember(model.rxns, dmemreaction)))=-100;
end


for i = 1:1:23
    disp(i)
    index = i
    index = num2str(index)
    disp(fileNames{i})
    filepath = strcat(folderPath,fileNames{i})
    data = readtable(filepath);
    expressionData = data.Var1;
    thre_ub = quantile(expressionData, 0.7)
    thre_lb = 0
    core_rxns = {'BIOMASS_maintenance', 'FACOAL40im', 'ACSm'};
    tol = 1e-6;
    cellModel = iMAT(model, expressionData,thre_lb, thre_ub, tol, core_rxns);
    sample_name = strsplit(fileNames{i},'_expression.csv')
    sample_name = strsplit(sample_name{1}, 'RxnScore_')
    sample_name = sample_name{2}
    save_gem_path = strcat('C:\Users\Alvis\Box\Boyu-Gut_microbiome\CyberGut\Case_Studies\Broad_Study\FAV_analysis_per_sample\model_imat_samples_0518_2024\Healthy_reconstructions\', sample_name, '.mat' );
    writeCbModel(cellModel, 'mat', save_gem_path);
end


% inflamed samples
folderPath = 'C:\Users\Alvis\Box\Boyu-Gut_microbiome\CyberGut\Case_Studies\Broad_Study\FAV_analysis_per_sample\Rxn_scores_0518_2024\Inflamed_RxnScore\'
files = dir(folderPath)
fileNames = {files(~[files.isdir]).name};
model = readCbModel('C:\Users\Alvis\Box\Boyu-Gut_microbiome\CyberGut\Draft_reconstructions\Draft_reconstruction_consensus\Upadate_05-08-2024\icolonoEpithelium_loopless.mat');
DMEM = table2cell(DMEMinputs);
model_EXCHANGE = strmatch('EX_', model.rxns);
model.lb(find(ismember(model.rxns,model.rxns(model_EXCHANGE))))=0;
model.ub(find(ismember(model.rxns,model.rxns(model_EXCHANGE))))=1000;

for i = 1:1:58
    dmemreaction = DMEM{i}
    model.lb(find(ismember(model.rxns, dmemreaction)))=-100;
end

for i = 1:1:13
    disp(i)
    index = i;
    index = num2str(index);
    disp(fileNames{i})
    filepath = strcat(folderPath,fileNames{i})
    data = readtable(filepath);
    expressionData = data.Var1;
    thre_ub = quantile(expressionData, 0.7)
    thre_lb = 0
    core_rxns = {'BIOMASS_maintenance', 'FACOAL40im', 'ACSm'};
    tol = 1e-6;
    cellModel = iMAT(model, expressionData,thre_lb, thre_ub, tol, core_rxns);
    sample_name = strsplit(fileNames{i},'_expression.csv')
    sample_name = strsplit(sample_name{1}, 'RxnScore_')
    sample_name = sample_name{2}
    save_gem_path = strcat('C:\Users\Alvis\Box\Boyu-Gut_microbiome\CyberGut\Case_Studies\Broad_Study\FAV_analysis_per_sample\model_imat_samples_0518_2024\Inflamed_reconstructions\', sample_name, '.mat' );
    writeCbModel(cellModel, 'mat', save_gem_path);
end



%% GEM per samples CD
% healthy samples
folderPath = 'C:\Users\Alvis\Box\Boyu-Gut_microbiome\CyberGut\Case_Studies\GSE164985\FAV_analysis_per_sample\Rxn_score_0518_2024\Healthy_Rxnscore\'

files = dir(folderPath)
fileNames = {files(~[files.isdir]).name};
model = readCbModel('C:\Users\Alvis\Box\Boyu-Gut_microbiome\CyberGut\Draft_reconstructions\Draft_reconstruction_consensus\Upadate_05-08-2024\icolonoEpithelium_loopless.mat');
DMEM = table2cell(DMEMinputs);
model_EXCHANGE = strmatch('EX_', model.rxns);
model.lb(find(ismember(model.rxns,model.rxns(model_EXCHANGE))))=0;
model.ub(find(ismember(model.rxns,model.rxns(model_EXCHANGE))))=1000;

for i = 1:1:58
    dmemreaction = DMEM{i}
    model.lb(find(ismember(model.rxns, dmemreaction)))=-100;
end



for i = 1:1:4
    disp(i)
    index = i;
    index = num2str(index);
    disp(fileNames{i})
    filepath = strcat(folderPath,fileNames{i})
    data = readtable(filepath);
    expressionData = data.Var1;
    thre_ub = quantile(expressionData, 0.7)
    thre_lb = 0
    core_rxns = {'BIOMASS_maintenance', 'FACOAL40im', 'ACSm'};
    tol = 1e-6;
    cellModel = iMAT(model, expressionData,thre_lb, thre_ub, tol, core_rxns);
    sample_name = strsplit(fileNames{i},'_expressionSet.csv')
    sample_name = strsplit(sample_name{1}, 'RxnScore_')
    sample_name = sample_name{2}
    save_gem_path = strcat('C:\Users\Alvis\Box\Boyu-Gut_microbiome\CyberGut\Case_Studies\GSE164985\FAV_analysis_per_sample\model_imat_0518_2024\Healthy\', sample_name, '.mat' );
    writeCbModel(cellModel, 'mat', save_gem_path);
end


% inflamed samples
folderPath = 'C:\Users\Alvis\Box\Boyu-Gut_microbiome\CyberGut\Case_Studies\GSE164985\FAV_analysis_per_sample\Rxn_score_0518_2024\Inflamed_Rxnscore\'
files = dir(folderPath)
fileNames = {files(~[files.isdir]).name};
model = readCbModel('C:\Users\Alvis\Box\Boyu-Gut_microbiome\CyberGut\Draft_reconstructions\Draft_reconstruction_consensus\Upadate_05-08-2024\icolonoEpithelium_loopless.mat');
DMEM = table2cell(DMEMinputs);
model_EXCHANGE = strmatch('EX_', model.rxns);
model.lb(find(ismember(model.rxns,model.rxns(model_EXCHANGE))))=0;
model.ub(find(ismember(model.rxns,model.rxns(model_EXCHANGE))))=1000;

for i = 1:1:58
    dmemreaction = DMEM{i}
    model.lb(find(ismember(model.rxns, dmemreaction)))=-100;
end



for i = 1:1:3
    disp(i)
    index = i;
    index = num2str(index);
    disp(fileNames{i})
    filepath = strcat(folderPath,fileNames{i})
    data = readtable(filepath);
    expressionData = data.Var1;
    thre_ub = quantile(expressionData, 0.7)
    thre_lb = 0
    core_rxns = {'BIOMASS_maintenance', 'FACOAL40im', 'ACSm'};
    tol = 1e-6;
    cellModel = iMAT(model, expressionData,thre_lb, thre_ub, tol, core_rxns);
    sample_name = strsplit(fileNames{i},'_expressionSet.csv')
    sample_name = strsplit(sample_name{1}, 'RxnScore_')
    sample_name = sample_name{2}
    save_gem_path = strcat('C:\Users\Alvis\Box\Boyu-Gut_microbiome\CyberGut\Case_Studies\GSE164985\FAV_analysis_per_sample\model_imat_0518_2024\Inflamed\', sample_name, '.mat' );
    writeCbModel(cellModel, 'mat', save_gem_path);
end


