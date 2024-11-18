load('workspaceThaueraAromatic.mat')
%lets load all the tables from Curation Thauera

oneRxnAdd=readcell("Curation Thauera.xlsx","Sheet",1);
oneRxnNo=readcell("Curation Thauera.xlsx","Sheet",2);
oneMetsAdd=readcell("Curation Thauera.xlsx","Sheet",3);
multRxnsAdd=readcell("Curation Thauera.xlsx","Sheet",4);
multRxnsNo=readcell("Curation Thauera.xlsx","Sheet",5);
metsMulAdd=readcell("Curation Thauera.xlsx","Sheet",6);
rxnsCarveMe=readcell("Curation Thauera.xlsx","Sheet",8);
MetsRxnsCM=readcell("Curation Thauera.xlsx","Sheet",9);
rxnsCarveMeNo=readcell("Curation Thauera.xlsx","Sheet",10);
rxnsCMIM=readcell("Curation Thauera.xlsx","Sheet",11);
RxnCMIMNo=readcell("Curation Thauera.xlsx","Sheet",12);

rxnsIn=struct('oneRxnAdd',{oneRxnAdd},'multRxnsAdd',{multRxnsAdd},...
    'rxnsCarveMe',{rxnsCarveMe},'rxnsCMIM',{rxnsCMIM});

rxnsOut=struct('oneRxnNo',{oneRxnNo},'multRxnsNo',{multRxnsNo},...
    'rxnsCarveMeNo',{rxnsCarveMeNo},'RxnCMIMNo',{RxnCMIMNo});

metsIn=struct('oneMetsAdd',{oneMetsAdd},'metsMulAdd',{metsMulAdd},...
    'MetsRxnsCM',{MetsRxnsCM});

%lets check which data in present in the model (should be rxnsIn and metsIn)
%and the reactions that should not be in the model (should be rxnsOut)
modelTest=modelTry6;
%present in model
fdNames = fieldnames(rxnsIn);
compTable = table();
for i = 1:length(fdNames)
    fieldName = fdNames{i};
    rxnInfo = rxnsIn.(fieldName);
    inModel = findRxnIDs(modelTest,rxnInfo(:,1));
    checkZero = find(inModel==0);
    rules = modelTest.grRules(inModel);
    rxnsIn.(fieldName)(:,7) = rules; 
    genes = findGenesFromRxns(modelTest,modelTest.rxns(inModel));
    containsFunc = @(str, substr) contains(str,substr);
    saveEqual = zeros(size(inModel,1),1); lenList = zeros(size(inModel,1),1);
    for j = 1:length(inModel)
        redList = genes{j,1};
        results = cellfun(@(x) arrayfun(@(s) containsFunc(x, s), redList), rules(j), 'UniformOutput', false);
        checkSize = isequal(size(results{1,1}==1,1),size(redList,1));
        saveEqual(j) = checkSize;
        lenList(j) = length(redList);
    end
    Tinfo = table(cellfun(@char,rxnInfo(:,1),"UniformOutput",false),lenList,saveEqual,'VariableNames',{'reactions','# of genes','match data'});
    compTable = vertcat(compTable,Tinfo);
end

%visualizing the data, all the grRules look correctly assigned and all the
%reactions are in the model
histogram(Tinfo.("match data"),2)
xticks([0, 1]);
%numerical validation
mistakes=isempty(find(Tinfo.("match data")==0, 1)); %empty!


%now lets check the reactions that should not be in the model
%since we are only interested in checking if the reactions are in the model
listNoRxns={oneRxnNo,multRxnsNo,rxnsCarveMeNo,RxnCMIMNo}';
rxnsNo = cellfun(@(z) z(:,1), listNoRxns,'UniformOutput', false);
listNo = vertcat(rxnsNo{:});
checkNo = findRxnIDs(modelTest,listNo); %some reactions are in the model!
%lets check if we have to delete these reactions, or if it's okay to keep
%them. This will happen in the gapfilling process
noIndex = checkNo(checkNo~=0);
rxnsToCheck = modelTest.rxns(noIndex); rxnsToCheck(:,2) = modelTest.grRules(noIndex);


%after that, lets verify the grRules and the remaining reactions to curate
%perform manual curation to check every single reaction and gene assigned

%take out the reactions from compTable column 1 and rxnsToCheck

rxSubstract = compTable.reactions;%rxns from compTable
rxSubstract = vertcat(rxSubstract,rxnsToCheck(:,1));

rxnsToCur = setdiff(modelTest.rxns,rxSubstract);
rxnIndex = findRxnIDs(modelTest,rxnsToCur);%indexes of these reactions
rxnsToCur(:,2) = modelTest.grRules(rxnIndex);

rxnsToCur = cell2table(rxnsToCur);
%lets create an excel file with the reactions and grRules of the remaining
%reactions

writetable(rxnsToCur,'curateThauera.xlsx','Sheet','DataToCurate')

%lets only keep metabolic reactions, removing exchange, demand and sink
%reactions


cellRxnsCur = table2cell(rxnsToCur);

eventsCur = cellfun(@(yz) yz(1:3), cellRxnsCur(:,1),'UniformOutput', false);
substringsToCheck = {'EX_', 'DM_', 'SK_'};
checkSubstring = @(str) any(contains(substringsToCheck, str));
pseudoRxns = ~cellfun(checkSubstring, eventsCur);

redCurThauera = cell2table(cellRxnsCur(pseudoRxns,:));
writetable(redCurThauera,'reducedCurThauera.xlsx','Sheet','redCuration')

%lets bring the information from the KEGGThauera, specifically we are
%interesed in the notBIGG file, which contains the ec numbers not found in
%bigg and we have to add manually.
%notBiGG contains only the EC numbers, we need to link the genes
%finalDict has the comparison of the NCBI and KEGG genes for Thauera
%ecGenes links the KEGG genes with the EC numbers 

%lets check if notBiGG and ecGenes contain unique EC numbers
uniqNotBiGG = unique(notBiGG); uniqEcGenes = unique(ecGenes(:,1));

for i=1:length(notBiGG)
    idxNot = find(strcmp(notBiGG{i,1},ecGenes(:,1)));
    uniqNotBiGG{i,2} = ecGenes(idxNot,2);
end

modNotBiGG = uniqNotBiGG(:,1);
for i=1:length(uniqNotBiGG)
    xGen = uniqNotBiGG{i,2};
    isMatch = find(cellfun(@(str) ismember(str, xGen), finalDict(:,1)));
    if ~isempty(isMatch)
        modNotBiGG{i,2} = finalDict(isMatch,2);
    else
        modNotBiGG{i,2} = [];
    end
end

%based on the EC numbers, lets get the rxn identifiers to link them to BiGG

webLink = 'https://rest.kegg.jp/get/ec:'; 
for i = 1: length(modNotBiGG)
    getRxnID = {};
    webData = webread(strcat(webLink,modNotBiGG{i,1}));
    rxnPos = strfind(webData,'ALL_REAC');
    if ~isempty(rxnPos)
        lines = strsplit(webData(rxnPos:end), '\n');
        getRxnID = strsplit(extractAfter(lines{1,1},'ALL_REAC    '));
        modNotBiGG{i,3} = getRxnID';
    else
        modNotBiGG{i,3} = [];
    end  
end

%lets get the ones that are not empty for the third column 

nonEmpty = find(~cellfun(@isempty,modNotBiGG(:,3)));
redNotBiGG = modNotBiGG(nonEmpty,:);



dataBIGG=readtable("bigg_models_reactions.xlsx","Sheet",1);
rxnKEGG =dataBIGG.bigg_id; rxnKEGG(:,2)=dataBIGG.database_links;

%lets look for the ones which have KEGG identifiers
txtKEGG = 'kegg.reaction/';
keggPos = strfind(rxnKEGG(:,2),txtKEGG);
keggIndex = find(~cellfun(@isempty,keggPos));
redKEGG = rxnKEGG(keggIndex,:);
redPos = keggPos(keggIndex);

for i=1:length(redKEGG)
    infoRxn = redPos{i,1};
    accumKEGG= {};
    for j=1:length(infoRxn)
        idxKEGG = erase(strtok(extractAfter(redKEGG{i,2}(infoRxn(j):end),txtKEGG)),';');
        accumKEGG{j,1} = idxKEGG;
    end
    redKEGG{i,3} = accumKEGG;
end

%lets see first any matching through a general intersect
rKEGG=redKEGG(:,3);keggFound=redNotBiGG(:,3);
keggBIGG = unique(vertcat(rKEGG{:})); keggNotBigg = unique(vertcat(keggFound{:}));
checkInt = intersect(keggBIGG,keggNotBigg);

%only 15 reactions were found in BIGG, lets find them
cont = 1;
posRXN = [];
for i = 1:numel(rKEGG)
    innerCellArray = rKEGG{i};
    innerIsMatch = cellfun(@(str) ismember(str, checkInt), innerCellArray, 'UniformOutput', false);
    numData = cell2mat(innerIsMatch);
    if ~isempty(find(numData==1, 1))
        i
        posRXN(cont,1) = i;
        cont = cont + 1;
    else
    end
end

foundRxns = redKEGG(posRXN);
connectedKEGG = redKEGG(posRXN,3);

webBigg='http://bigg.ucsd.edu/api/v2/universal/reactions/';

for i=1:length(foundRxns)
    rxnJson = webread(strcat(webBigg,foundRxns{i,1}));
    biggForm = strrep(rxnJson.reaction_string,'&#8652;','<=>');
    foundRxns{i,2} = biggForm;
end

%lets link the reactions to the genes
 for i = 1:length(connectedKEGG)
     innerKEGG = connectedKEGG{i,1};
     innerMatch = cellfun(@(str) ismember(str, innerKEGG), redNotBiGG(:,3), 'UniformOutput', false);
     sumElem = cellfun(@(innerCell) sum(innerCell), innerMatch);
     posElem = find(sumElem==1); infoElem = redNotBiGG(posElem,2);
     foundRxns{i,3}=unique(vertcat(infoElem{:}));
 end

%lets check if any of these reactions are already in the model
checkKEGG = findRxnIDs(modelTest,foundRxns(:,1));
rxnsNotModel = foundRxns(checkKEGG==0,:);
%rxnsNotModel were manually reviewed and classified in three groups in
%column E: value 0, reaction is not in the model and must be added with the
%reaction pointed in column D; value 1 is already in the model with a
%different reaction ID (reaction in model is listed in column D); and value
%2 which means non-equivalent reaction was found and there is no need to
%consider it 


curNotKEGG=readtable("reactionsToCheck.xlsx","Sheet",2);
curNotKEGG.Properties.VariableNames = {'reaction ID','reaction formula','genes','equiv rxn','value'};
%only add the ones with 0
filt0 = find(curNotKEGG.value==0);

rxnsAddKEGG = rxnsNotModel(filt0,:);
rxnsAddKEGG(:,4) = curNotKEGG.("equiv rxn")(filt0);

%now, lets work with the reactions that we couldn't find in BiGG using nor
%the EC number nor the KEGG ID

diffSet = setdiff(erase(keggNotBigg(3:end),';'),keggBIGG); %removing the noise of the 2 first rows

%lets map back and get the EC numbers we could not find in BIGG
%we will need redNotBiGG
cont = 1;
reviewKEGG={};
for i=1:length(redNotBiGG)
    keggInfo = redNotBiGG{i,3};
    matchKEGG = cellfun(@(str) ismember(str, keggInfo), diffSet, 'UniformOutput', false);
    getSum = sum(cell2mat(matchKEGG));
    if getSum>0
        reviewKEGG(cont,:) = redNotBiGG(i,:);
        cont = cont + 1;
    end
end

KEGGCheck = cell2table(reviewKEGG,"VariableNames",{'EC numbers','Genes Info','kegg IDs'});

%lets go back to the KEGG API 
for i = 1: length(reviewKEGG)
    webData = webread(strcat(webLink,reviewKEGG{i,1}));
    subPos = strfind(webData,'SUBSTRATE');
    prodPos = strfind(webData,'PRODUCT');
    comPos = strfind(webData,'COMMENT');
    lineSub = strsplit(webData(subPos:prodPos-1), '\n');
    lineSub = erase(lineSub,'SUBSTRATE  ');
    lineSub = erase(lineSub,' ');
    lineProd = strsplit(webData(prodPos:comPos-1), '\n');
    lineProd = erase(lineProd,'PRODUCT ');
    lineProd = erase(lineProd,' ');
    reviewKEGG{i,4} = lineSub';
    reviewKEGG{i,5} = lineProd';
 
end


%lets check the ones that contain both substrates and products
%%%MANUALLY REVIEWED: COULD BE IMPROVED USING CODE

reviewKEGG([51; 56; 59; 151; 152; 153; 154; 156; 181; 186],:) = [];

%lets get the CPD of the metabolites and check if we can find them in BIGG


metaBIGG=readtable("bigg_models_metabolites.xlsx","Sheet",1);

biggMets = metaBIGG.bigg_id; biggMets(:,2) = metaBIGG.database_links;
clear keggPos, clear keggIndex
txtComp = 'kegg.compound/';
keggPos = strfind(biggMets(:,2),txtComp);
keggIndex = find(~cellfun(@isempty,keggPos));
keggComp = biggMets(keggIndex,:);
redMets= keggPos(keggIndex);
 
for i=1:length(keggComp)
    infoMet = redMets{i,1};
    accumKEGG= {};
    for j=1:length(infoMet)
        idxKEGG = erase(strtok(extractAfter(keggComp{i,2}(infoMet(j):end),txtComp)),';');
        accumKEGG{j,1} = idxKEGG;
    end
    keggComp{i,3} = accumKEGG;
end

%lets extract from reviewKEGG the CPDs for substrate and products

detMets = reviewKEGG(:,1:2);
for i=1:length(reviewKEGG)

    modSub = reviewKEGG{i,4}(~cellfun(@isempty,reviewKEGG{i,4}));
    modProd = reviewKEGG{i,5}(~cellfun(@isempty,reviewKEGG{i,5}));
    extsub = erase(cellfun(@(str) strtok(str, '['), modSub,'UniformOutput',false),';');
    extProd = erase(cellfun(@(str) strtok(str, '['), modProd,'UniformOutput',false),';');
    cpdSub = erase(erase(erase(cellfun(@(str) extractAfter(str,'['),modSub,'UniformOutput',false),';'),']'),'CPD:');
    cpdProd = erase(erase(erase(cellfun(@(str) extractAfter(str,'['),modProd,'UniformOutput',false),';'),']'),'CPD:');
    detMets{i,3} = extsub; detMets{i,4} = extProd; detMets{i,5} = cpdSub; detMets{i,6} = cpdProd;
    detMets{i,7} = unique(vertcat(detMets{i,3},detMets{i,4})); detMets{i,7}=detMets{i,7}(~cellfun(@isempty,detMets{i,7}));     
    detMets{i,8} = unique(vertcat(detMets{i,5},detMets{i,6})); detMets{i,8}=detMets{i,8}(~cellfun(@isempty,detMets{i,8}));
end
lenCPD = cellfun(@length,detMets(:,8));
lenNames = cellfun(@length,detMets(:,7));
tableMets = table(detMets(:,1),detMets(:,2),lenNames,lenCPD,lenCPD./lenNames*100,...
    'VariableNames',{'EC Number','genes','# of metabolites','# KEGG CPDs','% of CPD found'});

%lets look for the intersection of the compounds from the ec numbers and
%the mets from bigg using the kegg identifiers

listCPD = vertcat(detMets{:,8});
interMets = intersect(listCPD,vertcat(keggComp{:,3}));
%lets get the BIGG IDs using the CPDs

matchKEGG2 = cellfun(@(str) ismember(str, interMets), keggComp(:,3), 'UniformOutput', false);
sums = arrayfun(@(sumVal) sum(cell2mat(sumVal)),matchKEGG2); 
posIdx = find(sums>0); posIdx2 = find(sums>1);
biggMetIDs = keggComp(posIdx,[1 3]); biggMetIDs2 = keggComp(posIdx2,1);
redIDs = cellfun(@(red) red(1:end-2),biggMetIDs(:,1),'UniformOutput',false);
uniqMets = unique(redIDs); lenUniqMets = length(uniqMets);

%lets compare detMets(:,8) with the interMets to chec

for i=1:length(detMets(:,8))
    keggInfo = detMets{i,8};
    matchKEGG = cellfun(@(str) ismember(str, keggInfo), interMets, 'UniformOutput', false);
    getSum = sum(cell2mat(matchKEGG));
    totSum(i,1) = getSum;   
    if getSum > 0
        detMets{i,9} = interMets(cell2mat(matchKEGG));
    else
        detMets{i,9} = [];
    end
    
end

tableMets.('mets matching BIGG') = totSum;
%lets work with the reduced bigg mets
%filter only the ones that are not empty in column 9
filt9 = find(~cellfun(@isempty,detMets(:,9)));
redDetMets = detMets(filt9,:);

for i=1:length(redDetMets)
    metsInv = redDetMets{i,9};
    cont = 1;
    modList = {};
    appendList = {};
    for j=1:length(biggMetIDs)
        xDict = biggMetIDs(j,:);
        modList = setdiff(xDict{1,2},modList);
        [intBigg,pos] = intersect(modList,metsInv);
        if ~isempty(intBigg)
            equivDict = intBigg; equivDict(:,2)=xDict(1);
            appendList = vertcat(appendList,equivDict);
            metsInv = setdiff(metsInv,appendList(:,1));
        end 
        if isempty(metsInv)
           break 
        end

    end
    appendList(:,2) = strcat(cellfun(@(red) red(1:end-2),appendList(:,2),'UniformOutput',false),'_c');
    redDetMets{i,10} = appendList;
end


%lets merge multiple models to check if we find all metabolites per
%reaction

load('iML1515.mat'),load('iAF692.mat'),load('iAF987.mat'),load('iJB785.mat')
load('iHN637.mat'),load('iYO844.mat'),load('iPC815.mat'),load('iSB619.mat')
load('iJN678.mat'),load('iJN746.mat'),load('iLJ478.mat'),load('iMM904.mat')
load('iRC1080.mat'),load('iIT341.mat'),load('iNF517.mat'),load('modelTest87Updated.mat')
load('Recon3D.mat'),load('iLB1027_lipid.mat')

models = {iML1515;iAF692;iAF987;iJB785;iHN637;iYO844;iPC815;iSB619;iJN678;iJN746;iLJ478;iMM904;iRC1080;...
    iIT341;iNF517;modelTest87;Recon3D;iLB1027_lipid};
listDesc = cellfun(@(field) getfield(field,'description'),models,'UniformOutput',false);
models = cellfun(@(struct,val) setfield(struct,'id',val),models, listDesc,'UniformOutput',false);

infoModels = cellfun(@(field) getfield(field,'rxns'),models,'UniformOutput',false);

%mergedModel=mergeModels(models);
resultCellArray = cellfun(@(list, str) [list'; repmat({str}, 1, numel(list))], infoModels, listDesc, 'UniformOutput', false)';
mergedRxns = horzcat(resultCellArray{:})';
[metsInv,stoichInv] = cellfun(@(mets,val) findMetsFromRxns(mets,val),models,infoModels,'UniformOutput',false);
fullMets = vertcat(metsInv{:});
[uniqRxns,uniqIdx] = unique(mergedRxns(:,1)); uniqRxns(:,2) = mergedRxns(uniqIdx,2);
uniqRxns(:,3) = fullMets(uniqIdx);

%lets remove the reactions with only 1 metabolite
%lets filter through length

getLength = cellfun(@length,uniqRxns(:,3));
idxMore1 = find(getLength>1); redRxns = uniqRxns(idxMore1,:);

%this for loop is gonna be very complex, we need to combine different
%situations. Probably this can be optimized with a different approach
for i=1:length(redDetMets)
    listMets = cellfun(@(red) red(1:end-2),redDetMets{i,10}(:,2),'UniformOutput',false);
    cont = 1;
    coinList = {};
    for j=1:length(uniqRxns)
        noComp = cellfun(@(red) red(1:end-2),uniqRxns{j,3},'UniformOutput',false);
        getInt = intersect(listMets,noComp);
        checkEq = isequaln(getInt,listMets);
        if checkEq
            coinList{cont,1} = uniqRxns{j,1};
            cont = cont + 1;              
        end
 
    end
    arrayCoinc{i,1} = coinList;
    redDetMets{i,11} = arrayCoinc{i,1};

end

redDetMets(:,11) = arrayCoinc;

tableRed = cell2table(redDetMets,'VariableNames',{'EC number','Genes involved','substrates','products','KEGG ID subs',...
    'KEGG ID prods','All metabolites','KEGG IDs identified','KEGG','mapping BIGG','rxns probably related'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Let's work with the model now
%check the quality of the simulations, the modifications in the BOF and
%gap-filling of metabolites involved in reactions with Thauera genes




