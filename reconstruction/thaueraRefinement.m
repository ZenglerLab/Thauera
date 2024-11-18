load('modelTryMod.mat')
[a,b] = exchangeSingleModel(modelTryMod);
model = modelTryMod;
initGrowth = optimizeCbModel(model);
fieldsRmv = {'metCharge';'metComps';'rxnFrom';'metFrom';'geneFrom';'rxnNotes';'rxnConfidenceScores'};
%remove non required fields
model = rmfield(model,fieldsRmv);
%reconstruct rev field
condition = model.lb<0;
model.rev = condition;
%lets close the multiple constraints above reactions 24 of the b vector
actvExch = b.A_Rxn; actvIdx = b.Index;
model.lb(actvIdx(25:end)) = 0;
model = changeRxnBounds(model,'EX_co2_e',0,'l');
model = changeRxnBounds(model,'EX_glc__D_e',-5,'l');
model = changeRxnBounds(model,'EX_o2_e',-5,'l');

[a1,b1] = exchangeSingleModel(model);
secGrowth = optimizeCbModel(model);

%lets check if we have repeated metabolites
draftModel1 = model;
uniqMets1 =unique(draftModel1.mets);
uniq1 = length(uniqMets1);

repMets1 = cell2mat(cellfun(@(mets) length(find(strcmp(draftModel1.mets,mets))),uniqMets1,'UniformOutput',false));
repIdx1 = cellfun(@(mets) find(strcmp(draftModel1.mets,mets)),uniqMets1,'UniformOutput',false);
idxRep1 = find(repMets1>1);
metsIdx1 = uniqMets1(idxRep1);
cleanIdx1 =repIdx1(idxRep1);
modelMod1 = draftModel1;
for i=1:length(metsIdx1)
    currIdx1 = cleanIdx1{i,1};
    getSt1 =find(draftModel1.S(currIdx1(end),:));
    val1 =full(draftModel1.S(currIdx1(end),getSt1));
    modelMod1.S(currIdx1(1),getSt1) = val1;
    modelMod1.S(currIdx1(end),getSt1) = 0;
   
end
%no duplicate metabolites!

%run mass balance to check reactions

[massImbalance, imBalancedMass, imBalancedCharge, imBalancedRxnBool, ~, ~, ~] = checkMassChargeBalance(model, 0);
% rmvMets1 = cell2mat(cellfun(@(idx) idx(end),cleanIdx1,'UniformOutput',false));
% modelMod1.S(rmvMets1,:)=[];
% modelMod1.metComps(rmvMets1) =[];
% modelMod1.b(rmvMets1) =[];
% modelMod1.mets(rmvMets1) =[];
% modelMod1.metNames(rmvMets1) =[];
% modelMod1.metFormulas(rmvMets1) =[];
% modelMod1.csense=modelMod1.csense(1:length(modelMod1.mets));
% modelMod1.metCharges = zeros(length(modelMod1.mets),1);

imbRed = find(~cellfun(@isempty,imBalancedMass));
massInv = imBalancedMass(imbRed);
rxnsInv = model.rxns(imbRed);
pseudoID = {'EX_';'DM_';'SK_';'BIOMASS'};
checkStrings = @(cellContent) any(contains(cellContent, pseudoID));
result = ~cellfun(checkStrings, rxnsInv);

rxnsUnb = rxnsInv(result); rxnsUnb(:,2) = massInv(result); rxnsUnb(:,3) = printRxnFormula(model,rxnsUnb(:,1));

[mets,~] = findMetsFromRxns(model,rxnsUnb(:,1));

getPos = cellfun(@(idx) findMetIDs(model,idx),mets,'UniformOutput',false);

getFormula = cellfun(@(formula) model.metFormulas(formula),getPos,'UniformOutput',false);

mergeInfo = cellfun(@(metForm,val) horzcat(metForm,val),mets,getFormula,'UniformOutput',false);

mergedList = horzcat(rxnsUnb,mergeInfo);

%lets bring the original data from bigg

webBase = 'http://bigg.ucsd.edu/api/v2/universal/reactions/';
%this code get problems at i=76, skip to 77
for i = 1 : length(rxnsUnb)
    webData = webread(strcat(webBase,rxnsUnb{i,1}));
    rxnFormula{i,1} = webData.reaction_string;
    rxnFormula{i,1} = strrep(rxnFormula{i,1},'&#8652;','=>');
    rxnFormula{i,2} = webData.metabolites;

end

mergedMod = horzcat(mergedList,rxnFormula);

%manual balancing 
rxnsFix = findRxnIDs(model,{'DABAAT';'DHNPA';'DINSK';'FOLR2'});
h_c = findMetIDs(model,'h_c');
model.S(h_c,rxnsFix(1)) = -1;
model.S(h_c,rxnsFix(2)) = 0;
dhglcn_c= findMetIDs(model,'2dhglcn_c');
model.metFormulas(dhglcn_c) = {'C6H9O7'};
actn__R_c= findMetIDs(model,'actn__R_c');
model.metFormulas(actn__R_c) = {'C4H8O2'};
diact_c= findMetIDs(model,'diact_c');
model.metFormulas(diact_c) = {'C4H6O2'};
alltt_c= findMetIDs(model,'alltt_c');
model.metFormulas(alltt_c) = {'C4H7N4O4'};
pmcoa_c= findMetIDs(model,'pmcoa_c');
model.metFormulas(pmcoa_c) = {'C28H41N7O19P3S'};
met2obut_c= findMetIDs(model,'4met2obut_c');
model.metFormulas(met2obut_c) = {'C5H8O3S'};
s_c= findMetIDs(model,'s_c');
model.metFormulas(s_c) = {'S'};
cmhm_c= findMetIDs(model,'5cmhm_c');
model.metFormulas(cmhm_c) = {'C8H5O7'};
cmhmsa_c= findMetIDs(model,'5cmhmsa_c');
model.metFormulas(cmhmsa_c) = {'C8H6O6'};
copre2_c = findMetIDs(model,'copre2_c');
model.metFormulas(copre2_c) = {'C42H38CoN4O16'};
copre3_c = findMetIDs(model,'copre3_c');
model.metFormulas(copre3_c) = {'C43H40CoN4O16'};
copre4_c = findMetIDs(model,'copre4_c');
model.metFormulas(copre4_c) = {'C44H43CoN4O16'};
copre5_c = findMetIDs(model,'copre5_c');
model.metFormulas(copre5_c) = {'C45H45CoN4O16'};
copre6_c = findMetIDs(model,'copre6_c');
model.metFormulas(copre6_c) = {'C44H45CoN4O16'};
codhpre6_c = findMetIDs(model,'codhpre6_c');
model.metFormulas(codhpre6_c) = {'C44H47CoN4O16'};
copre8_c = findMetIDs(model,'copre8_c');
model.metFormulas(copre8_c) = {'C45H52N4O14Co'};
fe3dhbzs_c = findMetIDs(model,'fe3dhbzs_c');
model.metFormulas(fe3dhbzs_c) = {'C10H10NO6Fe'};
fldox_c = findMetIDs(model,'fldox_c');
model.metFormulas(fldox_c) = {'X'};
fldox_c = findMetIDs(model,'fldrd_c');
model.metFormulas(fldox_c) = {'XH2'};
frmd_c = findMetIDs(model,'frmd_c');
model.metFormulas(frmd_c) = {'CH3NO'};
udpgalur_c = findMetIDs(model,'udpgalur_c');
model.metFormulas(udpgalur_c) = {'C15H22N2O18P2'};
ohed_c = findMetIDs(model,'2ohed_c');
model.metFormulas(ohed_c) = {'C7H6O5'};
dhhed_c = findMetIDs(model,'24dhhed_c');
model.metFormulas(dhhed_c) = {'C7H8O6'};
pprdn_c = findMetIDs(model,'pprdn_c');
model.metFormulas(pprdn_c) = {'C5H9N'};
rbl__D_c = findMetIDs(model,'rbl__D_c');
model.metFormulas(rbl__D_c) = {'C5H10O5'};
rbt_c = findMetIDs(model,'rbt_c');
model.metFormulas(rbt_c) = {'C5H12O5'};
sertrna_sec_c = findMetIDs(model,'sertrna_sec_c');
model.metFormulas(sertrna_sec_c) = {'C3H6NO2R'};
trnasecys_c = findMetIDs(model,'trnasecys_c');
model.metFormulas(trnasecys_c) = {'R'};
sucbz_c = findMetIDs(model,'sucbz_c');
model.metFormulas(sucbz_c) = {'C11H8O5'};
hba_c = findMetIDs(model,'4hba_c');
model.metFormulas(hba_c) = {'C7H8O2'};
tre_c = findMetIDs(model,'tre_c');
model.metFormulas(tre_c) = {'C12H22O11'};
indpyr_c = findMetIDs(model,'indpyr_c');
model.metFormulas(indpyr_c) = {'C11H8NO3'};
hethmpp_c = findMetIDs(model,'hethmpp_c');
model.metFormulas(hethmpp_c) = {'C14H20N4O8P2S'};
mdr1p_c = findMetIDs(model,'5mdr1p_c');
model.metFormulas(mdr1p_c) = {'C6H11O7PS'};
acmama_c = findMetIDs(model,'acmama_c');
model.metFormulas(acmama_c) = {'C14H23N2O9'};
ch15deccoa_c = findMetIDs(model,'ch15deccoa_c');
thfglu_c = findMetIDs(model,'thfglu_c');
model.metFormulas(thfglu_c) = {'C24H27N8O9'};
hh24dd_c = findMetIDs(model,'2hh24dd_c');
model.metFormulas(hh24dd_c) = {'C7H6O5'};
model.S(h_c,rxnsFix(3)) = 1;
model.S(h_c,rxnsFix(4)) = 1;
mdru1p_c = findMetIDs(model,'5mdru1p_c');
model.metFormulas(mdru1p_c) = {'C6H11O7PS'};
rnam_c = findMetIDs(model,'rnam_c');
model.metFormulas(rnam_c) = {'C11H15N2O5'};
hphaccoa_c = findMetIDs(model,'hphaccoa_c');
model.metFormulas(hphaccoa_c) = {'C29H38N7O18P3S'};

[massImbalance2, imBalancedMass2, imBalancedCharge2, imBalancedRxnBool2, ~, ~, ~] = checkMassChargeBalance(model, 0);


imbRed2 = find(~cellfun(@isempty,imBalancedMass2));
massInv2 = imBalancedMass2(imbRed2);
rxnsInv2 = model.rxns(imbRed2);
pseudoID = {'EX_';'DM_';'SK_';'BIOMASS'};
checkStrings2 = @(cellContent) any(contains(cellContent, pseudoID));
result2 = ~cellfun(checkStrings2, rxnsInv2);

rxnsUnb2 = rxnsInv2(result2); rxnsUnb2(:,2) = massInv2(result2); rxnsUnb2(:,3) = printRxnFormula(model,rxnsUnb2(:,1));

[mets2,~] = findMetsFromRxns(model,rxnsUnb2(:,1));

getPos2 = cellfun(@(idx) findMetIDs(model,idx),mets2,'UniformOutput',false);

getFormula2 = cellfun(@(formula) model.metFormulas(formula),getPos2,'UniformOutput',false);

mergeInfo2 = cellfun(@(metForm,val) horzcat(metForm,val),mets2,getFormula2,'UniformOutput',false);

mergedList2 = horzcat(rxnsUnb2,mergeInfo2);

%lets try to delete these reactions and check the impact in the 

growthModel = optimizeCbModel(model);
model2 = model; modelStock = model;
cont=1; cont2 =1;
for i = 1: length(mergedList2)
    model2 =removeRxns(model2,mergedList2{i,1});
    calcGrowth = optimizeCbModel(model2);
    if calcGrowth.f>0.01
        disp('reaction removed')
        modelStock = model2;
        rxnsRem{cont,1} = mergedList2{i,1};
        cont = cont + 1;
    else
        disp('reaction not removed')        
        model2=modelStock;
        rxnsKept{cont2,1} = mergedList2{i,1};
        cont2 = cont2 + 1;
    end
end

%lets look for complex metabolites that could be messing with the model
%lets check first these compounds: fmnRD_c, fmn_c, fdxo_42_c, fdxr_42_c,
%fdxo_2_2_c, fdxrd_c, mqn8_c, ficytb_c, focytc_c, fdxox_c,focytC_c, ficytC_c.
%Be careful with the Denitrification process!!

metsCheck = {'fmnRD_c', 'fdxo_42_c', 'fdxr_42_c','fdxo_2_2_c',...
    'fdxrd_c', 'ficytb_c', 'focytc_c', 'fdxox_c','focytC_c', 'ficytC_c'}';


revMets = findRxnsFromMets(model2,metsCheck);
revMets(:,2) = printRxnFormula(model2,revMets(:,1));
revMets(:,3) = model2.grRules(findRxnIDs(model2,revMets(:,1)));

%lets change the fmnRD_c for fmnh2_c, check if these reactions are already
%in the model.

model2 = changeRxnBounds(model2,{'EX_o2_e'},0,'u');
[a2,b2] = exchangeSingleModel(model2);
model2Opt = optimizeCbModel(model2);
%something weird is happening with oxygen
%lets keep track
%o2_e,o2_p,o2_c

%o2_e
o2e = findRxnsFromMets(model2,{'o2_e'});
o2e(:,2) = printRxnFormula(model2,o2e(:,1));
o2e(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,o2e(:,1))));

%o2_p
o2p = findRxnsFromMets(model2,{'o2_p'});
o2p(:,2) = printRxnFormula(model2,o2p(:,1));
o2p(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,o2p(:,1))));


%o2_c
o2c = findRxnsFromMets(model2,{'o2_c'});
o2c(:,2) = printRxnFormula(model2,o2c(:,1));
o2c(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,o2c(:,1))));

model2 = changeRxnBounds(model2,{'O2t'},0,'b');
optimizeCbModel(model2)
[a3,b3] = exchangeSingleModel(model2);

%lets check ammonium, glucose and phosphate fluxes
%glc__D_e
glc__D_e = findRxnsFromMets(model2,{'glc__D_e'});
glc__D_e(:,2) = printRxnFormula(model2,glc__D_e(:,1));
glc__D_e(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,glc__D_e(:,1))));

%glc__D_p
glc__D_p = findRxnsFromMets(model2,{'glc__D_p'});
glc__D_p(:,2) = printRxnFormula(model2,glc__D_p(:,1));
glc__D_p(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,glc__D_p(:,1))));


%glc__D_c
glc__D_c = findRxnsFromMets(model2,{'glc__D_c'});
glc__D_c(:,2) = printRxnFormula(model2,glc__D_c(:,1));
glc__D_c(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,glc__D_c(:,1))));

model2 = changeRxnBounds(model2,{'GLCtex_copy1';'GLCtex_copy2'},0,'b');
model2Opt=optimizeCbModel(model2);
[a4,b4] = exchangeSingleModel(model2);

%%%
%lets check again the oxygen to measure the new changes 

%o2_e
o2e_2 = findRxnsFromMets(model2,{'o2_e'});
o2e_2(:,2) = printRxnFormula(model2,o2e_2(:,1));
o2e_2(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,o2e_2(:,1))));

%o2_p
o2p_2 = findRxnsFromMets(model2,{'o2_p'});
o2p_2(:,2) = printRxnFormula(model2,o2p_2(:,1));
o2p_2(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,o2p_2(:,1))));


%o2_c
o2c_2 = findRxnsFromMets(model2,{'o2_c'});
o2c_2(:,2) = printRxnFormula(model2,o2c_2(:,1));
o2c_2(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,o2c_2(:,1))));

%glc__D_e
glc__D_e2 = findRxnsFromMets(model2,{'glc__D_e'});
glc__D_e2(:,2) = printRxnFormula(model2,glc__D_e2(:,1));
glc__D_e2(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,glc__D_e2(:,1))));

%glc__D_p
glc__D_p2 = findRxnsFromMets(model2,{'glc__D_p'});
glc__D_p2(:,2) = printRxnFormula(model2,glc__D_p2(:,1));
glc__D_p2(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,glc__D_p2(:,1))));


%glc__D_c
glc__D_c2 = findRxnsFromMets(model2,{'glc__D_c'});
glc__D_c2(:,2) = printRxnFormula(model2,glc__D_c2(:,1));
glc__D_c2(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,glc__D_c2(:,1))));


%malthx_c
malthx_c = findRxnsFromMets(model2,{'malthx_c'});
malthx_c(:,2) = printRxnFormula(model2,malthx_c(:,1));
malthx_c(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,malthx_c(:,1))));

%maltttr_c
maltttr_c = findRxnsFromMets(model2,{'maltttr_c'});
maltttr_c(:,2) = printRxnFormula(model2,maltttr_c(:,1));
maltttr_c(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,maltttr_c(:,1))));

%change the lower bounds of these reactions 'MLTG1','MLTG5','MLTG3','MLTG4'

rxnsGlc ={'MLTG1','MLTG5','MLTG3','MLTG4'};
model2 = changeRxnBounds(model2,rxnsGlc,0,'l');
model2Opt=optimizeCbModel(model2);
[a5,b5] = exchangeSingleModel(model2);

% change the lower bounds of these reactions 'DESAT16','DESAT14','DESAT18'

rxnsDesat ={'DESAT16','DESAT14','DESAT18'};
model2 = changeRxnBounds(model2,rxnsDesat,0,'l');
model2Opt=optimizeCbModel(model2);
[a6,b6] = exchangeSingleModel(model2);

model2=changeRxnBounds(model2,{'DHORDi'},0,'l');


%fe2_e
fe2_e = findRxnsFromMets(model2,{'fe2_e'});
fe2_e(:,2) = printRxnFormula(model2,fe2_e(:,1));
fe2_e(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,fe2_e(:,1))));

%fe2_p
fe2_p = findRxnsFromMets(model2,{'fe2_p'});
fe2_p(:,2) = printRxnFormula(model2,fe2_p(:,1));
fe2_p(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,fe2_p(:,1))));
%fe2_c
fe2_c = findRxnsFromMets(model2,{'fe2_c'});
fe2_c(:,2) = printRxnFormula(model2,fe2_c(:,1));
fe2_c(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,fe2_c(:,1))));

model2 = changeRxnBounds(model2,'FE2Gabcpp',0,'l');
model2Opt=optimizeCbModel(model2);
[a7,b7] = exchangeSingleModel(model2);

%fe3_e
fe3_e = findRxnsFromMets(model2,{'fe3_e'});
fe3_e(:,2) = printRxnFormula(model2,fe3_e(:,1));
fe3_e(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,fe3_e(:,1))));

%fe3_p
fe3_p = findRxnsFromMets(model2,{'fe3_p'});
fe3_p(:,2) = printRxnFormula(model2,fe3_p(:,1));
fe3_p(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,fe3_p(:,1))));
%fe3_c
fe3_c = findRxnsFromMets(model2,{'fe3_c'});
fe3_c(:,2) = printRxnFormula(model2,fe3_c(:,1));
fe3_c(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,fe3_c(:,1))));

%the KEY is in the fe3+!!!!

rxnsFe3 = {'FEENTERexs','FE3HOXexs','FECRMexs','ARBTNexs','FEOXAMexs','CPGNexs','FE3DCITexs','FERIRDe','FE3Gabcpp'};
model2 = changeRxnBounds(model2,rxnsFe3,0,'b');
model2Opt=optimizeCbModel(model2);
[a8,b8] = exchangeSingleModel(model2);

%nh4_e
nh4_e = findRxnsFromMets(model2,{'nh4_e'});
nh4_e(:,2) = printRxnFormula(model2,nh4_e(:,1));
nh4_e(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,nh4_e(:,1))));

%nh4_p
nh4_p = findRxnsFromMets(model2,{'nh4_p'});
nh4_p(:,2) = printRxnFormula(model2,nh4_p(:,1));
nh4_p(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,nh4_p(:,1))));
%nh4_c
nh4_c = findRxnsFromMets(model2,{'nh4_c'});
nh4_c(:,2) = printRxnFormula(model2,nh4_c(:,1));
nh4_c(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,nh4_c(:,1))));
model2 = changeRxnBounds(model2,{'NH4t','NH4t4pp'},0,'b');
model2Opt=optimizeCbModel(model2);
[a9,b9] = exchangeSingleModel(model2);


model2 = changeRxnBounds(model2,{'CYSTGL','GLUDxi','HAMR','ADD','CDGS','ATPHs','CSND','CTPS1','GTPHs','SERD_D'},0,'l');
model2 = changeRxnBounds(model2,{'NH3c'},0,'l');
model2Opt=optimizeCbModel(model2);
[a10,b10] = exchangeSingleModel(model2);


%h2o_e
h2o_e = findRxnsFromMets(model2,{'h2o_e'});
h2o_e(:,2) = printRxnFormula(model2,h2o_e(:,1));
h2o_e(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,h2o_e(:,1))));

%h2o_p
h2o_p = findRxnsFromMets(model2,{'h2o_p'});
h2o_p(:,2) = printRxnFormula(model2,h2o_p(:,1));
h2o_p(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,h2o_p(:,1))));

%h2o_c
h2o_c = findRxnsFromMets(model2,{'h2o_c'});
h2o_c(:,2) = printRxnFormula(model2,h2o_c(:,1));
h2o_c(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,h2o_c(:,1))));

rxnsH2Oe = {'H2Ot','TMAOR2e','TMAOR1e','GGTAe2'};
rxnsH2Op = {'CLPNH120pp','CLPNH140pp','CLPNH141pp','CLPNH160pp','CLPNH161pp','CLPNH180pp','CLPNH181pp'};
model2 = changeRxnBounds(model2,rxnsH2Oe,0,'b');
model2 = changeRxnBounds(model2,rxnsH2Op,0,'l');
model2Opt=optimizeCbModel(model2);
[a11,b11] = exchangeSingleModel(model2);

rxnsFix = {'ACOAD20','ACOADH2','AHMMPS_1','DHPDO','HPA3MO','FDMOtau','LYSMO','MTOLDOX','MXMO','PTOLDOX','PXMO','PYDXNO'};

rxnsClose = {'FDMO2_1','FDMO3_1','FDMO5_1','FDMO6_1','FDMO_1','O2t'};
model2 = changeRxnBounds(model2,rxnsClose,0,'b');
model2 = changeRxnBounds(model2,rxnsFix,0,'l');
model2Opt=optimizeCbModel(model2);
[a12,b12] = exchangeSingleModel(model2);


%so4_e
so4_e = findRxnsFromMets(model2,{'so4_e'});
so4_e(:,2) = printRxnFormula(model2,so4_e(:,1));
so4_e(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,so4_e(:,1))));

%so4_p
so4_p = findRxnsFromMets(model2,{'so4_p'});
so4_p(:,2) = printRxnFormula(model2,so4_p(:,1));
so4_p(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,so4_p(:,1))));

%so4_c
so4_c = findRxnsFromMets(model2,{'so4_c'});
so4_c(:,2) = printRxnFormula(model2,so4_c(:,1));
so4_c(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,so4_c(:,1))));

rxnsSO4 = {'SO4t2','SULabc'};
model2 = changeRxnBounds(model2,rxnsSO4,0,'b');
model2Opt=optimizeCbModel(model2);
[a13,b13] = exchangeSingleModel(model2);

%pi_e
pi_e = findRxnsFromMets(model2,{'pi_e'});
pi_e(:,2) = printRxnFormula(model2,pi_e(:,1));
pi_e(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,pi_e(:,1))));

%pi_p
pi_p = findRxnsFromMets(model2,{'pi_p'});
pi_p(:,2) = printRxnFormula(model2,pi_p(:,1));
pi_p(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,pi_p(:,1))));

%pi_c
pi_c = findRxnsFromMets(model2,{'pi_c'});
pi_c(:,2) = printRxnFormula(model2,pi_c(:,1));
pi_c(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,pi_c(:,1))));

rxnsPie = {'PIabc','PIt2r','PIt7'};
rxnsPip = {'3NTD2pp','3NTD9pp','3NTD4pp','3NTD7pp','2PGt6pp','3PGt6pp','PAPA181pp'};
model2 = changeRxnBounds(model2,rxnsPie,0,'b');
model2 = changeRxnBounds(model2,rxnsPip,0,'l');
model2Opt=optimizeCbModel(model2);
[a14,b14] = exchangeSingleModel(model2);


%h_e
h_e = findRxnsFromMets(model2,{'h_e'});
h_e(:,2) = printRxnFormula(model2,h_e(:,1));
h_e(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,h_e(:,1))));
%h_p
h_p = findRxnsFromMets(model2,{'h_p'});
h_p(:,2) = printRxnFormula(model2,h_p(:,1));
h_p(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,h_p(:,1))));
%h_c
h_c = findRxnsFromMets(model2,{'h_c'});
h_c(:,2) = printRxnFormula(model2,h_c(:,1));
h_c(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,h_c(:,1))));

rxnsHe = {'CADVt','ALAt2r','CD2t4','CRO4t3','D_LACt2','GLYt2r','HCINNMt2r','Kt2r','Kt3r','LYSt3r','NADHDH',...
   'NAt3_1','NO2t2r','NO3t2','NTR3B','PROt2r','RNF','SPMDt3','THD2','URATEt_1','URAt2','XANt2','ZN2t4',...
   'CHLt2','CITt12','DALAt2r','ADNt2pp_copy2','URIt2pp_copy2','GLCtex_copy1','GLCtex_copy2','THMDt2pp_copy2'};
model2 = changeRxnBounds(model2,rxnsHe,0,'b');
rxnsHc = {};

model2Opt=optimizeCbModel(model2);
[a15,b15] = exchangeSingleModel(model2);

%asn__L_e
asn__L_e = findRxnsFromMets(model2,{'asn__L_e'});
asn__L_e(:,2) = printRxnFormula(model2,asn__L_e(:,1));
asn__L_e(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,asn__L_e(:,1))));

%asn__L_p
asn__L_p = findRxnsFromMets(model2,{'asn__L_p'});
asn__L_p(:,2) = printRxnFormula(model2,asn__L_p(:,1));
asn__L_p(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,asn__L_p(:,1))));

%asn__L_c
asn__L_c = findRxnsFromMets(model2,{'asn__L_c'});
asn__L_c(:,2) = printRxnFormula(model2,asn__L_c(:,1));
asn__L_c(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,asn__L_c(:,1))));




%atp_e
atp_e = findRxnsFromMets(model2,{'atp_e'});
atp_e(:,2) = printRxnFormula(model2,atp_e(:,1));
atp_e(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,atp_e(:,1))));
%atp_p
atp_p = findRxnsFromMets(model2,{'atp_p'});
atp_p(:,2) = printRxnFormula(model2,atp_p(:,1));
atp_p(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,atp_p(:,1))));
%atp_c
atp_c = findRxnsFromMets(model2,{'atp_c'});
atp_c(:,2) = printRxnFormula(model2,atp_c(:,1));
atp_c(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,atp_c(:,1))));

rxnsAtp = {'PFK','PFK_2','PFK_3','PPDK','ADNCYC','AACOAT','FACOAL180','FACOAL1821','FACOAL141',...
   'FACOAL160','FACOAL80','THFGLUS','GLUCYS','PPCOAC','MCCC','ADADir','AIRC1','ANHMK','BCOALIG2',...
   'CA2abcpp','CRO4abcpp','Cut1','DGNSK','DINSK','FACOAL140','FACOAL161','FACOAL181','FACOAL40It2pp',...
   'FACOAL40t2pp','FACOAL50It2pp','GSNK','HEX4','INSCR','NADS2','PACCOAL2','PACCOAL3','PPCK',...
   'PYDXK','RBK','RNMK','SALCHS4FEabcpp','SALCHS4abcpp','SO3abcpp','TMDK1','URIK1'};

rxnsRvAtp = {'ALAabc','BUTSabc','ETHSabc','GLUabc','HEXSabc','ILEabc','ISTNTabc','LEUabc','MALTabc','METSRabc',...
    'METabc','PTRCabc','SUCCabc','THRabc','TSULabc','VALabc'};
model2 = changeRxnBounds(model2,rxnsRvAtp,0,'b');
model2 = changeRxnBounds(model2,rxnsAtp,0,'l');
model2Opt=optimizeCbModel(model2)
[a16,b16] = exchangeSingleModel(model2);

%definitely there are plenty of problems with the reversibility of several
%reactions after adding ecoli reactions
%model2 after BOF

BOF = find(model.c);

checkRxns = model2.rxns(BOF+1:end);
checkRxns(:,2) = printRxnFormula(model2,checkRxns(:,1));
checkRxns(:,3) = num2cell(model2.lb((BOF+1:end)));
checkRxns(:,4) = num2cell(model2.ub((BOF+1:end)));
checkRxns(:,5) = num2cell(model2Opt.x((BOF+1:end)));

%model2Stock=model2; STOCK model of model2 if something happens
%'SHSL4r' check this one later, CYS this one too, check NADHNQR, check 'NAD_H2'
%be careful with PFK_ppi, FUMAH, ADAPAT, ASPTA6, C120SN, C140SN, C160SN,
%DAPDA, DCMPDA, DKMPPD2, DPPS, ENLIPAtex, FDMO1, FMNRy, G3PD4, HOPNTAL3, ILEDHr
%IOR2b, IOR3b, IORb, MECDPDH4, MMSDHir, NADPHQR, NO3R3pp, OMLT, POR3b, SALCHS4FER3
%SALCHS4FEtonex, SALCHS4tex, TETtex, TETDHpp3, VOR1b, VOR3b
%check AM6SAD, 4M2OPLOXRD, 3HOPCD, ACACT6r_1, BTS3r, CYSDSF, FRD5, OCBT_1, ORLT
%PYDXO, SCYSSL, SSCOARy, 
rxnsFixLb = {'AHSERL2','AKGDa','AKGDb','OIVD1r','PDHcr','SHSL2r','SUCD1','UHBZ1t_pp','Q23DO',...
    'ASR2','PC20M','CODSCL5BMT','MTAP','VPAMTr','MAN1PT','PPCSCT','CRTNh','OPAH','CODSCL5DA',...
    'FUMAC','4OD2','CPH4S','HOPNTAL2','OP4ENH2','CDPGLC46DH','BCOALIG',...
   'COCHL','ALCD4y','3MCAT23DOX','4MCAT23DOX','CAT23DOX','3HCINNMH',...
   'HMSD2','HMSD','MMSAD2','3M2OPLOXRD','GTMLT','PRUK','FACOAE100',...
   'FACOAE140','FACOAE160','FACOAE180','FACOAE120','ACPpds','MALT',...
   'SUCR','ARGN','ALLTAHr','ADPRDP','UDPGALPpp','UDPGPpp','UACGALPpp',...
   'UACGAMPpp','ACYP','HMSH','HMSH2','PEPCK_re','UDPGLDC','PTHPS','GUACYC',...
   'MMM2','BACCL','PPCOAC','2HH24DDH1','3HBCOAHL','3HOXPACt2pp','3MBZALDH','3MBZDH',...
   '3MBt2pp','3MBt4pp','4HBALDt2pp','4HOXPACMOF','4HOXPACt2pp','4MBZALDH','4MBZDH',...
   '4OD','ACM6PH','ACOAD2','ACONCtupp','ACPS','ADMDC','ADNK3','AGMT','AGPAT140',...
   'AGPAT141','ALDD1','ALDD20x','ALDD3','ALDD31','ALDD31_1','ALDD4','ALDD6','ALR3',...
   'AMMQT8_2','AMPMS','APENTAMAH','APH120','APH140','APH141','APH160','APH161','APH180',...
   'APH181','ARGN_1','BUTt4pp','CDP4D6DGLCRx','CHOLSH','CLBtpp','CMHMI','CRO4t3pp',...
   'CYSS','CYTDH','DH3MCHCDH','DH4MCHCDH','DHEDAA','DHNAOT','DHNCOAT','DHNPA_1','DPHAPC100',...
   'DPHAPC120','DPHAPC121','DPHAPC140','DPHAPC141','DPHAPC60','DPHAPC80','ECOAH12','ECOAH9ir',...
   'FACOAE161','FOMETRi','G3PD6','G3PD7','GDPMNH','GLXO1','HDCAt2pp','HMPK2','HMPK3','HMPK4',...
   'HMSH3','HYPOE','INS2D','INSK','K2L4Atex','KAS16','LPLIPAL2E141','MACCOAT','MCITD',...
   'MCSNAH','MDRPD','METOX1s','METOX2s','MGSA','MI1PS','MLDEP1pp','MLDEP2pp','NADPHQR3',...
   'NMNHYD','NO3t2pp','NORZpp','NTP7','OAAt2_2pp','OGLT','OPTCCL','ORLT','OXOAEL','PAPA120',...
   'PAPA140','PAPA141','PAPA160','PAPA180','PETNT161pp','PETNT181pp','PPAt2pp','PPPGO2',...
   'PSSA140','PSSA141','PYRDC','PYRt4pp','P_XYLtpp','RAFHpp','SBP','SCYSSL','SPTc','UDPDPS',...
   'UM4PCP','VNTDM','XPPT'};

rxnsElim = {'G6PI','SER_AL','URAt2pp_copy2','CCP','Q23DO','PHYTES','PFK_ppi','G1PCTYT','AMAA','FUMAH',...
   'COA1819ZD9DS','2PGLYCt6','ACALDt','ADAPAT','ADNK4','ALAt4','ASPTA6', 'C120SN', 'C140SN',...
   'C160SN','CDPABEQS','CMLDC','CO2t','COAt','COBALTt5','CPK1','CPPPGOAN2', 'DAPDA','DASYN_EC',...
   'DCMPDA','DKMPPD2','DPPS','ENLIPAtex','ETHAt','FDMO1','FMNRy','FRUt3','G3PD4','GLCAASE3',...
   'GLCURT','GLYOX_1','H2CO3_NAt_syn','HACOADr','HOPNTAL3','ILEDHr','IOR2b','IOR3b','IORb',...
   'LALDO','LPLIPA1','MECDPDH4','MGt5','MMSDHir','NADPHQR','NARK','NH3c','NO3R3pp','NO3t7',...
   'NOt','OACT','OALT','OAO5t3ex','OAO5t3pp','OMLT','PASYN_EC','PGPP_EC','PGSA_EC',...
   'PLIPA1','PNTOt2','POR3b','PYDXNtr','SALCHS4FER3','SALCHS4FEtonex','SALCHS4tex',...
   'SERt4', 'SULAabc','St','TETDHpp3','TETtex','THPAT','THRt4','URAt','UREAt','VOR1b',...
   'VOR3b'};

model2 = changeRxnBounds(model2,rxnsFixLb,0,'l');
model2 = changeRxnBounds(model2,rxnsElim,0,'b');

model2Opt=optimizeCbModel(model2,'max','one',false)
[a17,b17] = exchangeSingleModel(model2);

posLb = find(model2.lb<-1000);
posUb = find(model2.ub>1000);

model2 = changeRxnBounds(model2,model2.rxns(posLb),-1000,'l');
model2 = changeRxnBounds(model2,model2.rxns(posUb),1000,'u');

activeRxns = model2.rxns(model2Opt.x~=0);
activeRxns(:,2) = printRxnFormula(model2,activeRxns(:,1));
activeRxns(:,3) = num2cell(model2.lb(model2Opt.x~=0));
activeRxns(:,4) = num2cell(model2.ub(model2Opt.x~=0));
activeRxns(:,5) = num2cell(model2Opt.x(model2Opt.x~=0));



[minFlux, maxFlux] = fluxVariability(model2, 'rxnNameList', {'EX_o2_e'})
[solution] = enumerateOptimalSolutions(model2);

model2 = changeRxnBounds(model2,{'EX_glc__D_e'},-3,'l');
model3Opt=optimizeCbModel(model2)
[a18,b18] = exchangeSingleModel(model2);
activeRxns2 = model2.rxns(model3Opt.x~=0);
activeRxns2(:,2) = printRxnFormula(model2,activeRxns2(:,1));
activeRxns2(:,3) = num2cell(model2.lb(model3Opt.x~=0));
activeRxns2(:,4) = num2cell(model2.ub(model3Opt.x~=0));
activeRxns2(:,5) = num2cell(model3Opt.x(model3Opt.x~=0));

model2 = changeRxnBounds(model2,{'EX_glc__D_e'},0,'b');
rxnsReview = {'SUCD4','NAD_H2'};
model2 = changeRxnBounds(model2,{rxnsReview{1,1},'BZt1pp'},0,'l');
model2 = changeRxnBounds(model2,rxnsReview(2),0,'b');
model2 = changeRxnBounds(model2,{'EX_co2_e'},-1000,'l');
model2Opt=optimizeCbModel(model2)
[a19,b19] = exchangeSingleModel(model2);
load('iML1515.mat')
[ec1,ec2] = exchangeSingleModel(iML1515);
rxnsToCur = {'CYSDSF','POR5','CYTDt2pp_copy2','NADHNQR','Natex','HSDx','HSDxi','SPT_syn',...
    'ARAB_Dt'};
model2 = changeRxnBounds(model2,rxnsToCur,0,'b');

%co2_e
co2_e = findRxnsFromMets(model2,{'co2_e'});
co2_e(:,2) = printRxnFormula(model2,co2_e(:,1));
co2_e(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,co2_e(:,1))));
%co2_p
co2_p = findRxnsFromMets(model2,{'co2_p'});
co2_p(:,2) = printRxnFormula(model2,co2_p(:,1));
co2_p(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,co2_p(:,1))));
%co2_c
co2_c = findRxnsFromMets(model2,{'co2_c'});
co2_c(:,2) = printRxnFormula(model2,co2_c(:,1));
co2_c(:,3) = num2cell(model2Opt.x(findRxnIDs(model2,co2_c(:,1))));

%lets check duplicate reactions 
uniqRxns = unique(model2.rxns);


%lets create an algorithm
[mets,~] = findMetsFromRxns(model2,model2.rxns);
cont = 1;
for i = 1 :length(model2.rxns)
    otherElements = mets([1:i-1, i+1:end]);
    rxnsList = model2.rxns([1:i-1, i+1:end]);
    checkEq = find(cell2mat(cellfun(@(equal) isequaln(mets{i,1},equal),otherElements,'UniformOutput',false)));
    if ~isempty(checkEq)
        rxnsSim{cont,1} = model2.rxns{i,1};
        rxnsSim{cont,2} = printRxnFormula(model2,model.rxns{i,1});
        rxnsSim{cont,3} = rxnsList(checkEq);
        rxnsSim{cont,4} = printRxnFormula(model2,rxnsList(checkEq));
        cont = cont +1 ;
    end
end

%lets check which reactions we can remove from the list, since we are only
%repeating some reactions 
rxnsDel = {'HACD5','GLYOX_1','HACD6i','HACD4i','HACD7i','HACD1i','HACD2i','HOPNTAL3',...
    'CYSS2','CYSabc2pp','ZN2abcpp','HACD3i','EX_pnto__e','URAt2pp_copy1','HACD8i',...
    'ACOAD4_1','ACACT5r_1','SERR','VALDHr','GAPDi_nadp','INSt3pp','KAT1', 'ACOAD1fr'};
%review 'INSt3pp', 'KAT1', 'ACOAD1fr'
model2 = changeRxnBounds(model2,rxnsDel,0,'b');


model3Opt=optimizeCbModel(model2)
[a20,b20] = exchangeSingleModel(model2);
activeRxns3 = model2.rxns(model3Opt.x~=0);
activeRxns3(:,2) = printRxnFormula(model2,activeRxns3(:,1));
activeRxns3(:,3) = num2cell(model2.lb(model3Opt.x~=0));
activeRxns3(:,4) = num2cell(model2.ub(model3Opt.x~=0));
activeRxns3(:,5) = num2cell(model3Opt.x(model3Opt.x~=0));


model2 = changeRxnBounds(model2,{'EX_glc__D_e'},-3,'l');
model3Opt=optimizeCbModel(model2)
[a21,b21] = exchangeSingleModel(model2);
activeRxns4 = model2.rxns(model3Opt.x~=0);
activeRxns4(:,2) = printRxnFormula(model2,activeRxns4(:,1));
activeRxns4(:,3) = num2cell(model2.lb(model3Opt.x~=0));
activeRxns4(:,4) = num2cell(model2.ub(model3Opt.x~=0));
activeRxns4(:,5) = num2cell(model3Opt.x(model3Opt.x~=0));


model3 = model2;
rxnsTotRmv = {rxnsDel',{'GLCtex_copy1','GLCtex_copy2'}',rxnsFe3',{'NH3c'},rxnsClose',...
    rxnsPie',rxnsHe',rxnsRvAtp',rxnsElim',rxnsToCur',rxnsDel'}';
rmvRxns = unique(vertcat(rxnsTotRmv{:}));
model3 =removeRxns(model3,rmvRxns);
model4Opt=optimizeCbModel(model3)
[a22,b22] = exchangeSingleModel(model3);


%lets check which metabolites are not being used in any reaction 
metS=size(model3.S,1);
cont = 1;
for i = 1 : metS
    getStoich = find(full(model3.S(i,:))~=0);
    if isempty(getStoich)
        metRmv{cont,1} = model3.mets{i,1};
        metRmv{cont,2} = {i};
        cont = cont +1;
    end
end

model3 = removeMetabolites(model3,metRmv(:,1));


%lets now remove the exchange reactions and metabolites linked which are
%not involved in the experimental data and the biolog plates,
%PM1
PM1 = readcell('Supplementary Material 2.xlsx','Sheet',1);
exchBiolog1 = PM1(2:80,3);
redExch1 = unique(exchBiolog1(contains(exchBiolog1,'_e')));
%PM2
PM2 = readcell('Supplementary Material 2.xlsx','Sheet',2);
exchBiolog2 = PM2(2:65,3);
redExch2 = unique(exchBiolog2(contains(exchBiolog2,'_e')));
%PM3
PM3 = readcell('Supplementary Material 2.xlsx','Sheet',3);
exchBiolog3 = PM3(2:90,2);
redExch3 = unique(exchBiolog3(contains(exchBiolog3,'_e')));

AAs = {'arg__L_e';'asn__L_e';'asp__L_e';'cys__L_e';'gln__L_e';'glu__L_e';...
    'gly_e';'his__L_e';'ile__L_e';'leu__L_e';'lys__L_e';'met__L_e';'ala__L_e';...
    'phe__L_e';'pro__L_e';'ser__L_e';'thr__L_e';'trp__L_e';'tyr__L_e';...
    'val__L_e'};

exchData = strcat('EX_',unique(vertcat(redExch1,redExch2,redExch3,AAs)));
exchModel = model3.rxns(contains(model3.rxns,'EX_'));

exchInt = intersect(exchModel,exchData);
exchMissing = setdiff(exchData,exchModel);
exchExtra = setdiff(exchModel,exchData);

exchCheck = exchExtra;
exchCheck(:,2) = printRxnFormula(model3,exchCheck(:,1));
exchCheck(:,3) = model3.metNames(findMetIDs(model3,erase(exchCheck(:,1),'EX_')));
%writecell(exchCheck,'ExchangeToCheck.xlsx')

%lets bring the filtered information 

filtExch = readcell('ExchangeToCheck.xlsx','Sheet',1);
filtExch2 = readtable('ExchangeToCheck.xlsx','Sheet',1);
filtVal = filtExch2.Var4;
metP = strcat(erase(erase(filtExch(:,1),'EX_'),'_e'),'_p');
rxnsP = cellfun(@(mets) findRxnsFromMets(model3,mets),metP,'UniformOutput',false);
metE = erase(filtExch(:,1),'EX_');
metE2 = metE(filtVal==1);
rxnsE = cellfun(@(mets) findRxnsFromMets(model3,mets),metE,'UniformOutput',false);
[check1,check2] = exchangeSingleModel(model3);
exchNot = erase(check2.A_Rxn,'EX_');
redMetE = setdiff(metE2,exchNot);

%for loop to delete reactions which are not killing the model
model4 = model3;
modelPiv = model4; cont = 1;
removedMets = {};
for i = 1: length(redMetE)
   model4 = removeRxns(model4,findRxnsFromMets(model4,redMetE{i,1}));
   growthVal = optimizeCbModel(model4);
   if growthVal.f > 0.001
       modelPiv = model4;
       removedMets{cont,1} = redMetE{i,1}; 
       cont = cont + 1;
   else
       model4 = modelPiv; 
   end
end

metS2=size(model4.S,1);
cont = 1;
metRmv = {};
for i = 1 : metS2
    getStoich = find(full(model4.S(i,:))~=0);
    if isempty(getStoich)
        metRmv{cont,1} = model4.mets{i,1};
        metRmv{cont,2} = {i};
        cont = cont +1;
    end
end

%Lets check the connectivity inside the model
%lets look for dead ends 

thaueraCur = readcell('curateThauera1.xlsx','Sheet',1);
redCur  = thaueraCur(2:end,:);
%lets check which reactions are still in the model
remainCur = findRxnIDs(model4,redCur(:,1));
idxRem = remainCur(remainCur~=0);
delRxns = redCur(remainCur==0);
remRxns = redCur(remainCur~=0);
remIdx2 = cell2mat(cellfun (@(str) strmatch(str,redCur(:,1),'exact'),remRxns,'UniformOutput',false));
remData = redCur(remIdx2,:);

modRules = cellfun(@ismissing,remData(:,3),'UniformOutput',false);
calcLen = cellfun(@length,modRules);
rulesWithEmpty = remData(:,3); idxMiss = find(calcLen==1);
rulesWithEmpty(idxMiss,1) = {''};
rulesRed = cellfun(@char,rulesWithEmpty,'UniformOutput',false);
checkEmpty = find(cellfun(@isempty,rulesRed));

%lets check which grRules suposed to have a gene (even an exogenous gene)
uncRules = remData(:,2);
modRules2 = cellfun(@ismissing,uncRules,'UniformOutput',false);
calcLen2 = cellfun(@length,modRules2);
idxMiss2 = find(calcLen2==1);
uncRules(idxMiss2,1) = {''};
rulesRed2 = cellfun(@char,uncRules,'UniformOutput',false);
checkEmpty2 = find(cellfun(@isempty,rulesRed2));
checkEcol = find(contains(uncRules,'b'));

intEmptyEcol = intersect(checkEmpty,checkEcol);

rulesMerge = horzcat(remData(intEmptyEcol,1),uncRules(intEmptyEcol),rulesRed(intEmptyEcol));

idtRxn = findRxnIDs(model4,rulesMerge(:,1));

%lets try to run a blast with the ecoli genes against the TMZ1T genes
eraseRules = erase(erase(erase(erase(rulesMerge(:,2),'('),')'),'or'),'and');
splitInfo = cellfun(@(genes) strsplit(genes)',eraseRules,'UniformOutput',false);
bGenes = unique(vertcat(splitInfo{:}));
webColi = 'http://bigg.ucsd.edu/api/v2/models/iML1515/genes/';
for i = 1 : length(bGenes)
    infoGenes = webread(strcat(webColi,bGenes{i,1}));
    redEcoli{i,1} = bGenes{i,1};
    redEcoli{i,2} = infoGenes.protein_sequence;
end

%fastawrite('redEcoli.fasta',redEcoli(:,1),redEcoli(:,2))

ecoliThauera = readcell('ecoliThauera.csv');

VarNames = {'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',...
    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'}';
tableRed = table(ecoliThauera(:,1),ecoliThauera(:,2),ecoliThauera(:,3),...
    ecoliThauera(:,4),ecoliThauera(:,11),ecoliThauera(:,12),'VariableNames',{VarNames{1},...
    VarNames{2},VarNames{3},VarNames{4},VarNames{11},VarNames{12}});

%eval filter
eval =1e-4;
cellRed = table2cell(tableRed);
filtEval = find(cell2mat(cellRed(:,5))<eval);
filtBlast = cellRed(filtEval,:);

%lets get the best hits

listChar=cellfun(@char,filtBlast(:,1),'UniformOutput',false);
[list,A,C]=unique(listChar);
clear bestHit
for i=1:length(list)
    pos=find(strcmp(listChar,list{i,1}));
    identVal=cell2mat(filtBlast(pos,3));
    proof=filtBlast(pos,1);
    [getMax,idxMax]=max(identVal);
    bestHit(i,:)=filtBlast(pos(idxMax),:);
    
end


rulesMerge2 = horzcat(rulesMerge,model4.rxnNames(findRxnIDs(model4,rulesMerge(:,1))));
rulesMerge2(:,5) = rulesMerge2(:,2);
for i = 1 : length(bestHit)
    rulesMerge2(:,5)= strrep(rulesMerge2(:,5),bestHit{i,1},bestHit{i,2});
end

%lets curate the remaining ecoli genes manually, but first let's see which
%ones are actually used in other reactions for _p mets

idxMerge2 = rulesMerge2(contains(rulesMerge2(:,5),'b'),:);
idxMerge2(:,6) = printRxnFormula(model4,idxMerge2(:,1));
%Lets bring the periplasm curation 
curationP = readcell('curationPeriplasm.xlsx','Sheet',1);
modRules3 = cellfun(@ismissing,curationP(:,7),'UniformOutput',false);
calcLen3 = cellfun(@length,modRules3);
idxMiss3 = find(calcLen3==1);
curationP(idxMiss3,7) = {''};
rulesRed3 = cellfun(@char,curationP(:,7),'UniformOutput',false);
nonEmpty = unique(rulesRed3(~cellfun(@isempty,rulesRed3)));

%check which of the nonEmpty list match with the exchange metabolites
%removed in the previous step
%metsE contains the exchangeMetabolites removed from the model through
%Exchange Reactions
metsMod = strrep(nonEmpty,'_p','_e');
intEP = intersect(metsMod,removedMets);
metsMod2 = strrep(intEP,'_e','_p');
rxnsP2 = cellfun(@(mets) findRxnsFromMets(model4,mets),metsMod2,'UniformOutput',false);
rxnsP2full = findRxnsFromMets(model4,metsMod2);



%for loop to delete reactions which are not killing the model
model5 = updateGenes(model4);
model5 = generateRules(model5);
modelPiv = model5; cont = 1;
removedRxns = {};
for i = 1: length(rxnsP2full)
   model5 = removeRxns(model5,rxnsP2full{i,1});
   growthVal = optimizeCbModel(model5);
   if growthVal.f > 0.001
       modelPiv = model5;
       removedRxns{cont,1} = rxnsP2full{i,1}; 
       cont = cont + 1;
   else
       model5 = modelPiv; 
   end
end

metS3=size(model5.S,1);
cont = 1;
metRmv = {};
for i = 1 : metS3
    getStoich = find(full(model5.S(i,:))~=0);
    if isempty(getStoich)
        metRmv{cont,1} = model5.mets{i,1};
        metRmv{cont,2} = {i};
        cont = cont +1;
    end
end

model5 = generateRules(model5);
model5 = updateGenes(model5);


%now lets update the grRules that we know have been corrected from file
%CurateThauera1

%lets check which reactions are still in the model
remainCur = findRxnIDs(model5,redCur(:,1));
idxRem = remainCur(remainCur~=0);
delRxns = redCur(remainCur==0);
remRxns = redCur(remainCur~=0);
remIdx2 = cell2mat(cellfun (@(str) strmatch(str,redCur(:,1),'exact'),remRxns,'UniformOutput',false));
remData = redCur(remIdx2,:);

modRules = cellfun(@ismissing,remData(:,3),'UniformOutput',false);
calcLen = cellfun(@length,modRules);
rulesWithEmpty = remData(:,3); idxMiss = find(calcLen==1);
rulesWithEmpty(idxMiss,1) = {''};
rulesRed = cellfun(@char,rulesWithEmpty,'UniformOutput',false);
checkEmpty = find(cellfun(@isempty,rulesRed));

curToAdd = horzcat(remRxns,rulesRed);

idxModel5 = findRxnIDs(model5,remRxns);
model5.grRules(idxModel5) = curToAdd(:,2);
model5 = buildRxnGeneMat(model5);
model5 = generateRules(model5);
model5 = updateGenes(model5);

%Finally, lets check dead ends
%check point
checkPoint5 = optimizeCbModel(model5);
[chp1,chp2] = exchangeSingleModel(model5);
DEmodel5 = detectDeadEnds(model5);
modelDE = model5.mets(DEmodel5);
rxnsDE = cellfun(@(mets) findRxnsFromMets(model5,mets),modelDE,'UniformOutput',false);
mergeDE = horzcat(modelDE,rxnsDE);
formDE = cellfun(@(rxns) printRxnFormula(model5,rxns),rxnsDE,'UniformOutput',false);
mergeDE(:,3) = formDE;
rulesDE = cellfun(@(rxns) model5.grRules(findRxnIDs(model5,rxns)),rxnsDE,'UniformOutput',false);
mergeDE(:,4) = rulesDE;



%lets get exchange reactions and test the production of each metabolite

allExch5 = model5.rxns(contains(model5.rxns,'EX_'));
%get Carbon metabolites
cont3 = 1;
for i = 1 : length(allExch5)
    model6 = changeRxnBounds(model5,allExch5{i,1},-2,'l');
    optGrowth = optimizeCbModel(model6);
    cont = 1; cont2 = 1; unfMets = {}; feasMets ={};
    i
    if optGrowth.f > 0.001
        for j = 1 : length(modelDE)
            model7 = addSinkReactions2(model6,modelDE{j,1});
            model7 = changeRxnBounds(model7,model7.rxns(end),0.01,'b');
            opt6 = optimizeCbModel(model7);
            if opt6.f>0.001
                feasMets{cont,1} = modelDE{j,1};
                cont = cont +1;
            else
                unfMets{cont2,1} = modelDE{j,1};
                cont2 = cont2 +1;
            end
        end
        feasData{cont3,1} = feasMets;
        unfData{cont3,1} = unfMets;
        cont3 = cont3 +1;
    end
end


%try adding all exchange
cont3 = 1; model6 = model5;
for i = 1 : length(allExch5)
    model6 = changeRxnBounds(model6,allExch5{i,1},-2,'l');
    optGrowth = optimizeCbModel(model6);
    cont = 1; cont2 = 1; unfMets2 = {}; feasMets2 ={};
    i
    if optGrowth.f > 0.001
        for j = 1 : length(modelDE)
            model7 = addSinkReactions2(model6,modelDE{j,1});
            model7 = changeRxnBounds(model7,model7.rxns(end),0.01,'b');
            opt6 = optimizeCbModel(model7);
            if opt6.f>0.001
                feasMets2{cont,1} = modelDE{j,1};
                cont = cont +1;
            else
                unfMets2{cont2,1} = modelDE{j,1};
                cont2 = cont2 +1;
            end
        end
        feasData2{cont3,1} = feasMets2;
        unfData2{cont3,1} = unfMets2;
        cont3 = cont3 +1;
    end
end


uniqFeasMets2 = unique(vertcat(feasData2{:}));
%uniqFeasMets2 are metabolites that can be sinthetyzed if exchange
%reactions are opened

%lets get the list of reactions which even opening the exchanges we don't
%get any growth

remDE = setdiff(modelDE,uniqFeasMets2);
rxnsDE2 = cellfun(@(mets) findRxnsFromMets(model5,mets),remDE,'UniformOutput',false);
mergeDE2 = horzcat(remDE,rxnsDE2);
formDE2 = cellfun(@(rxns) printRxnFormula(model5,rxns),rxnsDE2,'UniformOutput',false);
mergeDE2(:,3) = formDE2;
rulesDE2 = cellfun(@(rxns) model5.grRules(findRxnIDs(model5,rxns)),rxnsDE2,'UniformOutput',false);
mergeDE2(:,4) = rulesDE2;

findEmpty = find(cell2mat(cellfun(@(genes) sum(contains(genes,'TMZ1T')),mergeDE2(:,4),'UniformOutput',false))==0);
redMergeDE = mergeDE2(findEmpty,:);
model6= model5;
modelPiv = model5; cont = 1;
removedRxns2 = {}; 
for i = 1: length(redMergeDE)
   model6 = removeRxns(model6,redMergeDE{i,2});
   growthVal = optimizeCbModel(model6);
   if growthVal.f > 0.001
       modelPiv = model6;
       removedRxns2{cont,1} = redMergeDE{i,2}; 
       cont = cont + 1;
   else
       model6 = modelPiv; 
   end
end

metS4=size(model6.S,1);
cont = 1;
metRmv = {};
for i = 1 : metS4
    getStoich = find(full(model6.S(i,:))~=0);
    if isempty(getStoich)
        metRmv{cont,1} = model6.mets{i,1};
        metRmv{cont,2} = {i};
        cont = cont +1;
    end
end

model6 = generateRules(model6);
model6 = updateGenes(model6);

%lets check each remaining deadend

[remDE2,idxDE] = setdiff(mergeDE2(:,1),redMergeDE(:,1));
mergeDE3 = mergeDE2(idxDE,:);
collapsedDE = vertcat(mergeDE3{:,2});
collapsedDE(:,2) = vertcat(mergeDE3{:,3});
collapsedDE(:,3) = vertcat(mergeDE3{:,4});
%writecell(collapsedDE,'deadEndsCheck.csv')

find(contains(model6.rxns,'23DOGULNt4pp'))

curDE = readtable('deadEndsCheck2.xlsx','Sheet',1);

[rxnDE, idxDE2] = unique(curDE.Var1);
rxnsCritDE = curDE.Var4(idxDE2);

rxnsToDel = rxnDE(rxnsCritDE==0);
%delete dead end reactions obtained from the excel file deadEndsCheck2.xlsx

model7= model6;
modelPiv = model6; cont = 1;
removedRxns2 = {}; 
for i = 1: length(rxnsToDel)
   model7 = removeRxns(model7,rxnsToDel{i,1});
   growthVal = optimizeCbModel(model7);
   if growthVal.f > 0.001
       modelPiv = model7;
       removedRxns2{cont,1} = rxnsToDel{i,1}; 
       cont = cont + 1;
   else
       model7 = modelPiv; 
   end
end

[test1,test2] = exchangeSingleModel(model7);
metS5=size(model7.S,1);
cont = 1;
metRmv = {};
for i = 1 : metS5
    getStoich = find(full(model7.S(i,:))~=0);
    if isempty(getStoich)
        metRmv{cont,1} = model7.mets{i,1};
        metRmv{cont,2} = {i};
        cont = cont +1;
    end
end

model7 = buildRxnGeneMat(model7);
model7 = generateRules(model7);
model7 = updateGenes(model7);
model7 = removeUnusedGenes(model7);
checkGenesB = model7.grRules(contains(model7.grRules,'b'));

%lets add the new reactions and metabolites from the excel file 

newRxns = readcell('deadEndsCheck2.xlsx','Sheet',2);

newMets = readcell('deadEndsCheck2.xlsx','Sheet',3);

modRules4 = cellfun(@ismissing,newRxns(:,4),'UniformOutput',false);
calcLen4 = cellfun(@length,modRules4);
rulesWithEmpty4 = newRxns(:,4); idxMiss4 = find(calcLen4==1);
rulesWithEmpty4(idxMiss4,1) = {''};
rulesRed4 = cellfun(@char,rulesWithEmpty4,'UniformOutput',false);
checkEmpty4 = find(cellfun(@isempty,rulesRed4));

newRxns(:,4) = rulesRed4;

model8 = addMultipleMetabolites(model7,newMets(:,1),'metNames',newMets(:,2),'metFormulas',newMets(:,3),'metCharges',newMets(:,4));

modelPiv = model8; cont = 1;
addedRxns = {}; 
for i = 1: length(newRxns)
   model8 = addReaction(model8,newRxns{i,1},'reactionName',newRxns{i,2},'reactionFormula',newRxns{i,3},'geneRule',newRxns{i,4});
   growthVal = optimizeCbModel(model8);
   if growthVal.f > 0.001
       modelPiv = model8;
       i
       addedRxns{cont,1} = newRxns{i,1}; 
       cont = cont + 1;
   else
       model8 = modelPiv; 
   end
end

[test1,test2] = exchangeSingleModel(model8);
closeBd = {'EX_3mb_e';'EX_udpacgal_e';'EX_adocbl_e';'EX_all__D_e';'EX_arbtn_fe3_e';...
    'EX_aso4_e';'EX_btoh_e';'EX_butso3_e';'EX_g3pc_e';'EX_chols_e';'EX_glyc2p_e';...
    'EX_ibt_e';'EX_mepn_e';'EX_pentso3_e'};

model8 = changeRxnBounds(model8,closeBd,0,'l');

model8 = buildRxnGeneMat(model8);
model8 = generateRules(model8);
model8 = updateGenes(model8);
model8 = removeUnusedGenes(model8);



%lets modify Biomass amino acids and trace elements
%https://link.springer.com/article/10.1007/s00253-014-5756-x

%copy files 
%Mineral requirements.xlsx, BiomassMineralsThauera.m, BiomassRhodo.m

model9=model8;
stoichTable = readtable('aminoacids thesis.xlsx','Sheet',3);
bofID = find(model9.c==1);
AAnamesCyto=strcat(stoichTable.aminoAcidID,'_c');
stoichAA = stoichTable.StoichiometricValues;
metsID=findMetIDs(model9,AAnamesCyto);

for i=1:length(metsID)
    oldValue(i,1)=model9.S(metsID(i),bofID);
    model9.S(metsID(i),bofID)=-stoichAA(i,1);
    newValue(i,1)=-stoichAA(i,1);
    comparison(i,1)=abs(100*(newValue(i,1)-oldValue(i,1))/oldValue(i,1));
    
end


mediaChange=mean(comparison);
AAchanged=AAnames(find((comparison==(max(comparison)))));

%Biomass minerals

[bofMets,bofStoich]=findMetsFromRxns(model9,model9.rxns(find(model9.c)));
bofMets=bofMets{1,1};bofStoich=full(bofStoich{1,1});
MWdata=computeMW2(model9,bofMets);


%load minerals
[stoichMinerals,totMinerals,~]=xlsread('Mineral requirements.xlsx','Global');
totMinerals=totMinerals(2:end,1:2);
totMinerals(:,3)=strrep(totMinerals(:,2),'_e','_c');
thaueraStoich=stoichMinerals(:,4);
idxAct=find(thaueraStoich~=0);
totMinerals=totMinerals(idxAct,1:3);
stoichRed=thaueraStoich(idxAct);
RedData=horzcat(totMinerals,num2cell(stoichRed));
metsCheck=[totMinerals(:,2); totMinerals(:,3)];
metsPT=findMetIDs(model9,metsCheck);
notInModel=metsCheck(find(metsPT==0));
bofPres=intersect(totMinerals(:,3),bofMets);
notInBof=setdiff(totMinerals(:,3),bofMets);

%add to the BOF and transporters
model10=model9;
na_pos=findMetIDs(model10,'na1_c');
model10.S(na_pos,find(model10.c))=-stoichRed(3,1);
%[e1,e2]=exchangeSingleModel(model10);
%model10=changeRxnBounds(model10,{'EX_o2_e'},-10,'l');
[e3,e4]=exchangeSingleModel(model10);
model10=addMetabolite(model10,'bo3_c','metName','Borate','metFormula','BO3');
model10=addMetabolite(model10,'bo3_e','metName','Borate','metFormula','BO3');
model10=addMetabolite(model10,'bo3_p','metName','Borate','metFormula','BO3');
model10=addExchangeRxn(model10,{'bo3_e'});
model10=addReaction(model10,'BORtex','reactionName','Borate transport via diffusion (extracellular to periplasm)','reactionFormula','bo3_e <=> bo3_p','subSystem','Transport');
model10=addReaction(model10,'BORtpp','reactionName','Borate transport via diffusion (periplasm to cytoplasm)','reactionFormula','bo3_p <=> bo3_c','subSystem','Transport');
bo3_pos=findMetIDs(model10,'bo3_c');
findMineral=findMetIDs(model10,totMinerals(:,3));
model10.S(findMineral,find(model10.c))=-stoichRed(:,1);
model10.S(bo3_pos,find(model10.c))=-stoichRed(11,1);
[e5,e6]=exchangeSingleModel(model10);
model10=changeRxnBounds(model10,{'EX_tungs_e';'EX_slnt_e';'EX_sel_e'},0,'l');
optimizeCbModel(model10)



%lets add the KEGG information curated by Camila
rxnsKEGG = readtable('nitrogenRxns.xlsx','Sheet',2);
metsKEGG = readtable('nitrogenRxns.xlsx','Sheet',3);

%lets remove repeated reactions that we are gonna add later for
%denitrification
%no2_p and no2_c
rxnsNo2c = findRxnsFromMets(model10,'no2_c');
rxnsNo2c(:,2) = printRxnFormula(model10,rxnsNo2c(:,1));
rxnsNo2p = findRxnsFromMets(model10,'no2_p');
rxnsNo2p(:,2) = printRxnFormula(model10,rxnsNo2p(:,1));
%no3_p no3_c
rxnsNo3c = findRxnsFromMets(model10,'no3_c');
rxnsNo3c(:,2) = printRxnFormula(model10,rxnsNo3c(:,1));
rxnsNo3p = findRxnsFromMets(model10,'no3_p');
rxnsNo3p(:,2) = printRxnFormula(model10,rxnsNo3p(:,1));

remNitrogen = {'NTRIR2x';'NO3R2pp';'NTRIR4pp';'NTRIR3pp';'NORZpp';'NO3R1pp';...
    'NO3R2bpp'};

%model10 = removeRxns(model10,remNitrogen);

%first lets add these new metabolites

metIDs = metsKEGG.Var1;
checkInMets = findMetIDs(model10,metIDs);
metsAdd = metIDs(checkInMets==0);
clear dataKEGG
dataKEGG(:,1) = metsKEGG.Var1; dataKEGG(:,2) = metsKEGG.Var2; 
dataKEGG(:,3) = metsKEGG.Var3; dataKEGG(:,4) = num2cell(metsKEGG.Var4); 
detailMets = dataKEGG(checkInMets==0,:);

%check which reactions are already in the model
clear keggInfo
checkInRxns = findRxnIDs(model10,rxnsKEGG.BiGGID);
keggInfo = rxnsKEGG.BiGGID;keggInfo(:,2) =rxnsKEGG.Name;
keggInfo(:,3) =rxnsKEGG.Formula; keggInfo(:,4) = rxnsKEGG.GrRule;
detailRxns = keggInfo(checkInRxns==0,:);
%exchange validated
%https://ami-journals.onlinelibrary.wiley.com/doi/full/10.1111/j.1462-2920.2004.00615.x?casa_token=dsTgvtPNllIAAAAA%3AZEO0BXNtrAuFIk5NsMDFhpbuF6w9VIBNBoWwVsn_pD9fKE8LTpodgneYNVASVIxAnWrgJYu06Rzktd7DlA
%citrate and succinate, acetate
uniqueMets = unique(detailMets(:,1));
counts = histcounts(categorical(detailMets(:,1)), uniqueMets);
repeated_elements = uniqueMets(counts > 1);
model10 = addMultipleMetabolites(model10,detailMets(:,1),'metNames',detailMets(:,2),'metFormulas',detailMets(:,3),'metCharges',detailMets(:,4));
model11 =model10;
modelPiv = model11; cont = 1;
addedRxns = {}; 
for i = 1: length(detailRxns)
   model11 = addReaction(model11,detailRxns{i,1},'reactionName',detailRxns{i,2},'reactionFormula',detailRxns{i,3},'geneRule',detailRxns{i,4});
   growthVal = optimizeCbModel(model11);
   if growthVal.f > 0.001
        growthVal.f
       modelPiv = model11;
       addedRxns{cont,1} = detailRxns{i,1}; 
       cont = cont + 1;
   else
       model11 = modelPiv; 
   end
end

addedRxns(:,2) = printRxnFormula(model11,addedRxns(:,1));

metS6=size(model11.S,1);
cont = 1;
metRmv = {};
for i = 1 : metS6
    getStoich = find(full(model11.S(i,:))~=0);
    if isempty(getStoich)
        metRmv{cont,1} = model11.mets{i,1};
        metRmv{cont,2} = {i};
        cont = cont +1;
    end
end
model11 = removeMetabolites(model11, metRmv(:,1));
model11 = changeRxnBounds(model11,'EX_co2_e',0,'l');
[e7,e8]=exchangeSingleModel(model11);
optimizeCbModel(model11)
model11 = changeRxnBounds(model11,'EX_n2_e',0,'l');


model12 = model11; 
model12 = changeRxnBounds(model12,{'EX_no3_e';'EX_no2_e'},-5,'l');
model12 = changeRxnBounds(model12,'EX_o2_e',-3,'l');
[e9,e10]=exchangeSingleModel(model12);
optData=optimizeCbModel(model12);

%no3_p no3_c
rxnsNo3c = findRxnsFromMets(model12,'no3_c');
rxnsNo3c(:,2) = printRxnFormula(model12,rxnsNo3c(:,1));
rxnsNo3c(:,3) =num2cell(optData.x(findRxnIDs(model12,rxnsNo3c(:,1))));
rxnsNo3p = findRxnsFromMets(model12,'no3_p');
rxnsNo3p(:,2) = printRxnFormula(model12,rxnsNo3p(:,1));
rxnsNo3p(:,3) =num2cell(optData.x(findRxnIDs(model12,rxnsNo3p(:,1))));


%q8
q8c = findRxnsFromMets(model12,'q8_c');
q8c(:,2) = printRxnFormula(model12,q8c(:,1));
q8c(:,3) =num2cell(optData.x(findRxnIDs(model12,q8c(:,1))));


%q8h2
q8h2c = findRxnsFromMets(model12,'q8h2_c');
q8h2c(:,2) = printRxnFormula(model12,q8h2c(:,1));
q8h2c(:,3) =num2cell(optData.x(findRxnIDs(model12,q8h2c(:,1))));


%no2_p and no2_c
rxnsNo2c = findRxnsFromMets(model12,'no2_c');
rxnsNo2c(:,2) = printRxnFormula(model12,rxnsNo2c(:,1));
rxnsNo2c(:,3) =num2cell(optData.x(findRxnIDs(model12,rxnsNo2c(:,1))));

rxnsNo2p = findRxnsFromMets(model12,'no2_p');
rxnsNo2p(:,2) = printRxnFormula(model12,rxnsNo2p(:,1));
rxnsNo2p(:,3) =num2cell(optData.x(findRxnIDs(model12,rxnsNo2p(:,1))));


%no_p and no_c
rxnsNoc = findRxnsFromMets(model12,'no_c');
rxnsNoc(:,2) = printRxnFormula(model12,rxnsNoc(:,1));
rxnsNoc(:,3) =num2cell(optData.x(findRxnIDs(model12,rxnsNoc(:,1))));

rxnsNop = findRxnsFromMets(model12,'no_p');
rxnsNop(:,2) = printRxnFormula(model12,rxnsNop(:,1));
rxnsNop(:,3) =num2cell(optData.x(findRxnIDs(model12,rxnsNop(:,1))));


%n2o_p and n2o_c
rxnsN2oc = findRxnsFromMets(model12,'n2o_c');
rxnsN2oc(:,2) = printRxnFormula(model12,rxnsN2oc(:,1));
rxnsN2oc(:,3) =num2cell(optData.x(findRxnIDs(model12,rxnsN2oc(:,1))));

rxnsN2op = findRxnsFromMets(model12,'n2o_p');
rxnsN2op(:,2) = printRxnFormula(model12,rxnsN2op(:,1));
rxnsN2op(:,3) =num2cell(optData.x(findRxnIDs(model12,rxnsN2op(:,1))));


%ficytc_c and focytc_c
rxnsfocytc_c = findRxnsFromMets(model12,'focytc_c');
rxnsfocytc_c(:,2) = printRxnFormula(model12,rxnsfocytc_c(:,1));
rxnsfocytc_c(:,3) =num2cell(optData.x(findRxnIDs(model12,rxnsfocytc_c(:,1))));

rxnsN2op = findRxnsFromMets(model12,'n2o_p');
rxnsN2op(:,2) = printRxnFormula(model12,rxnsN2op(:,1));
rxnsN2op(:,3) =num2cell(optData.x(findRxnIDs(model12,rxnsN2op(:,1))));

newRxnsAdd ={'CYTBMQOR2pp_2';'NO3R3pp_2';'BCR6pp'};
newRxnsAdd(:,2) = {'Cytochrome b Menaquinone oxidoreductase';...
   'Nitrate reductase (Cytochrome c) periplasmic';'BenzoylCoA reductase (Ech associated, periplasm) 6 protons';...
   };
newRxnsAdd(:,3) = {'mql8_c + 2.0 ficytc_c => mqn8_c + 2.0 h_p + 2.0 focytc_c';...
    '4.0 h_c + no3_c + 2.0 focytc_c => h2o_c + no2_c + 2.0 h_p + 2.0 ficytc_c';...
    '6.0 h_p + benzcoa_c + fdxrd_c <=> 4.0 h_c + ch15deccoa_c + fdxox_c'};

model13 =model12;
modelPiv = model13; cont = 1;
addedRxns = {}; 
for i = 1: size(newRxnsAdd,1)
   model13 = addReaction(model13,newRxnsAdd{i,1},'reactionName',newRxnsAdd{i,2},'reactionFormula',newRxnsAdd{i,3});
   growthVal = optimizeCbModel(model13);
   if growthVal.f > 0.001
        growthVal.f
       modelPiv = model13;
       addedRxns{cont,1} = newRxnsAdd{i,1}; 
       cont = cont + 1;
   else
       model13 = modelPiv; 
   end
end
[b1,b2] = exchangeSingleModel(model13);
addedRxns(:,2) = printRxnFormula(model13,addedRxns(:,1));
%model13 = changeRxnBounds(model13,'EX_o2_e',0,'l');
optData = optimizeCbModel(model13);
%no3_p no3_c
rxnsNo3c = findRxnsFromMets(model13,'no3_c');
rxnsNo3c(:,2) = printRxnFormula(model13,rxnsNo3c(:,1));
rxnsNo3c(:,3) =num2cell(optData.x(findRxnIDs(model13,rxnsNo3c(:,1))));
rxnsNo3p = findRxnsFromMets(model13,'no3_p');
rxnsNo3p(:,2) = printRxnFormula(model13,rxnsNo3p(:,1));
rxnsNo3p(:,3) =num2cell(optData.x(findRxnIDs(model13,rxnsNo3p(:,1))));


%q8
q8c = findRxnsFromMets(model13,'q8_c');
q8c(:,2) = printRxnFormula(model13,q8c(:,1));
q8c(:,3) =num2cell(optData.x(findRxnIDs(model13,q8c(:,1))));


%q8h2
q8h2c = findRxnsFromMets(model13,'q8h2_c');
q8h2c(:,2) = printRxnFormula(model13,q8h2c(:,1));
q8h2c(:,3) =num2cell(optData.x(findRxnIDs(model13,q8h2c(:,1))));


%no2_p and no2_c
rxnsNo2c = findRxnsFromMets(model13,'no2_c');
rxnsNo2c(:,2) = printRxnFormula(model13,rxnsNo2c(:,1));
rxnsNo2c(:,3) =num2cell(optData.x(findRxnIDs(model13,rxnsNo2c(:,1))));

rxnsNo2p = findRxnsFromMets(model13,'no2_p');
rxnsNo2p(:,2) = printRxnFormula(model13,rxnsNo2p(:,1));
rxnsNo2p(:,3) =num2cell(optData.x(findRxnIDs(model13,rxnsNo2p(:,1))));


%no_p and no_c
rxnsNoc = findRxnsFromMets(model13,'no_c');
rxnsNoc(:,2) = printRxnFormula(model13,rxnsNoc(:,1));
rxnsNoc(:,3) =num2cell(optData.x(findRxnIDs(model13,rxnsNoc(:,1))));

rxnsNop = findRxnsFromMets(model13,'no_p');
rxnsNop(:,2) = printRxnFormula(model13,rxnsNop(:,1));
rxnsNop(:,3) =num2cell(optData.x(findRxnIDs(model13,rxnsNop(:,1))));


%n2o_p and n2o_c
rxnsN2oc = findRxnsFromMets(model13,'n2o_c');
rxnsN2oc(:,2) = printRxnFormula(model13,rxnsN2oc(:,1));
rxnsN2oc(:,3) =num2cell(optData.x(findRxnIDs(model13,rxnsN2oc(:,1))));

rxnsN2op = findRxnsFromMets(model13,'n2o_p');
rxnsN2op(:,2) = printRxnFormula(model13,rxnsN2op(:,1));
rxnsN2op(:,3) =num2cell(optData.x(findRxnIDs(model13,rxnsN2op(:,1))));


%ficytc_c and focytc_c
rxnsfocytc_c = findRxnsFromMets(model13,'focytc_c');
rxnsfocytc_c(:,2) = printRxnFormula(model13,rxnsfocytc_c(:,1));
rxnsfocytc_c(:,3) =num2cell(optData.x(findRxnIDs(model13,rxnsfocytc_c(:,1))));

rxnsficytc_c = findRxnsFromMets(model13,'ficytc_c');
rxnsficytc_c(:,2) = printRxnFormula(model13,rxnsficytc_c(:,1));
rxnsficytc_c(:,3) =num2cell(optData.x(findRxnIDs(model13,rxnsficytc_c(:,1))));


%mql8_c
rxnsmql8_c = findRxnsFromMets(model13,'mql8_c');
rxnsmql8_c(:,2) = printRxnFormula(model13,rxnsmql8_c(:,1));
rxnsmql8_c(:,3) =num2cell(optData.x(findRxnIDs(model13,rxnsmql8_c(:,1))));


model13 = changeRxnBounds(model13,'EX_o2_e',0,'l');
model13 = changeRxnBounds(model13,'EX_co2_e',-10,'l');
optimizeCbModel(model13)
%lets activate multiple exchange reactions and see which one works
%except oxygen
exchTest = b1.A_Rxn;
findO2 = find(contains(exchTest,'EX_o2_e'));
exchTest2=exchTest([1:findO2-1 findO2+1:end],1);
%!curl -O 'http://bigg.ucsd.edu/static/models/iJN1463.mat'
model14 = model13;
for i = 1: size(exchTest2,1)
   model14 = changeRxnBounds(model14,exchTest2{i,1},-3,'l');
   growthVal = optimizeCbModel(model14);
   if growthVal.f == 0
        growthVal.f

   else
       i
       break
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

load('iAF987.mat')
modelTemp = iAF987;

rxnsGeo = setdiff(modelTemp.rxns,model13.rxns);
rxnsGeo(:,2) = modelTemp.rxnNames(findRxnIDs(modelTemp,rxnsGeo(:,1)));
rxnsGeo(:,3) = printRxnFormula(modelTemp,rxnsGeo(:,1));
model14 = model13;
modelPiv = model13; cont = 1;
addedRxns = {}; 
model14 = changeRxnBounds(model14,'EX_nh4_e',-5,'l');
for i = 1: size(rxnsGeo,1)
   model14 = addReaction(model14,rxnsGeo{i,1},'reactionName',rxnsGeo{i,2},'reactionFormula',rxnsGeo{i,3});
   growthVal = optimizeCbModel(model14);
   if growthVal.f == 0
        growthVal.f
       modelPiv = model14;
       addedRxns{cont,1} = rxnsGeo{i,1}; 
       cont = cont + 1;
       
   else
       i
       model14 = modelPiv; 
       break
   end
end

[exch1,exch2] = exchangeSingleModel(model14);

%check weird mets 
met1611 = model13.mets{1611,1};
rxns1611 = findRxnsFromMets(model13,met1611);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%exchange reactions
%modelStock= model14;
fluxGrow = optimizeCbModel(model14);
getAct = find(fluxGrow.x~=0);
rxnsAct = model14.rxns(getAct); 
rxnsAct(:,2) = printRxnFormula(model14,rxnsAct(:,1));
rxnsAct(:,3) = num2cell(fluxGrow.x(getAct));


%%%%%%%%%%%%%%%%%%%
%model13
modelAnr = changeRxnBounds(model13,'EX_o2_e',-0.2,'l');          
[an1,an2] = exchangeSingleModel(modelAnr);
growthAnr = optimizeCbModel(modelAnr);

%check Active Anoxic rxns
getAct2 = find(growthAnr.x~=0);
rxnsAct2 = modelAnr.rxns(getAct2); 
rxnsAct2(:,2) = printRxnFormula(modelAnr,rxnsAct2(:,1));
rxnsAct2(:,3) = num2cell(growthAnr.x(getAct2));
 
%check Active Aerobic rxns
modelAer = model13;
modelAer = changeRxnBounds(modelAer,'EX_no3_e',0,'l');
[aer1,aer2] = exchangeSingleModel(modelAer);
growthAer = optimizeCbModel(modelAer);
getAct3 = find(growthAer.x~=0);
rxnsAct3 = modelAer.rxns(getAct3); 
rxnsAct3(:,2) = printRxnFormula(modelAer,rxnsAct3(:,1));
rxnsAct3(:,3) = num2cell(growthAer.x(getAct3));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%identified from the anerobicRxns.m script
fixAnr = {'PDX5PO2'};
fixAnr(:,2) = {'Pyridoxine 5-phosphate oxidase (anaerboic)'};
fixAnr(:,3) = {'nad_c + pdx5p_c => h_c + nadh_c + pydx5p_c'};

modelAnr2 =  changeRxnBounds(modelAnr,'EX_o2_e',0,'l');
modelAnr2 =  changeRxnBounds(modelAnr2,'EX_o2_e',1000,'u');
checkAnr = optimizeCbModel(modelAnr2)

addedRxns3={}; cont = 1;
modelPiv = modelAnr2;
for i = 1: size(fixAnr,1)
   modelAnr2 = addReaction(modelAnr2,fixAnr{i,1},'reactionName',fixAnr{i,2},'reactionFormula',fixAnr{i,3});
   growthVal = optimizeCbModel(modelAnr2)
   if growthVal.f == 0
        growthVal.f
       modelPiv = modelAnr2;
       addedRxns3{cont,1} = rxnsGeo{i,1}; 
       cont = cont + 1;
       
   else
       i
       break
   end
end

[anr3,anr4] = exchangeSingleModel(modelAnr2);

%%%
%MODEL FINALLY WORKING!!!!1
%modelAnaerobic version ======> modelAnr2

%lets 
load('iDT1278.mat')
%lets trace back the PHB metabolites and reactions
checkPHB = find(contains(iDT1278.subSystems,'PHB'));
PHBHint = iDT1278.rxns(checkPHB);
[metsPHB, ~] = findMetsFromRxns(iDT1278,PHBHint);
PHBHint(:,2) = printRxnFormula(iDT1278,PHBHint(:,1));
PHBHint(:,3) = metsPHB; 

phbAzoto = findRxnsFromMets(iDT1278, 'PHB_c');
phbAzoto2 = findRxnsFromMets(iDT1278,'phbg_c');
phbAzoto3 = findRxnsFromMets(iDT1278,'3hbcoa__R_c');
phbAzoto4 = findRxnsFromMets(iDT1278,'aacoa_c');
phbAzoto4(:,2) = printRxnFormula(iDT1278,phbAzoto4(:,1));

rxnsToCheck = unique(vertcat(PHBHint(:,1),phbAzoto,phbAzoto2,phbAzoto3,phbAzoto4(:,1)));
rxnsToCheck(:,2) = iDT1278.rxnNames(findRxnIDs(iDT1278,rxnsToCheck(:,1)));
rxnsToCheck(:,3) = printRxnFormula(iDT1278,rxnsToCheck(:,1));
[phbRel, ~] = findMetsFromRxns(iDT1278,rxnsToCheck(:,1));
rxnsToCheck(:,4) = phbRel;

%check in Thauera Model 

checkThauera = findRxnIDs(modelAnr2,rxnsToCheck(:,1));
rxnsFound = rxnsToCheck(checkThauera~=0,:);
rxnsMissing = rxnsToCheck(checkThauera==0,:);

%check which reactions 
wholeMissing = rxnsMissing(:,4);
metsToVer = unique(vertcat(wholeMissing{:}));
verfMets = findMetIDs(modelAnr2,metsToVer);
metsNotFound = metsToVer(verfMets==0);
metsFound = metsToVer(verfMets~=0);

%can we produce all the mets found in both models?
modelAer2 = modelAnr2;
modelAer2 = changeRxnBounds(modelAer2,'EX_o2_e',-0.5,'l');
growthAer = optimizeCbModel(modelAer2);
[ar1,ar2] = exchangeSingleModel(modelAer2);
for i = 1 : length(metsFound)
    modelAerCheck = addSinkReactions(modelAer2,metsFound{i,1});
    modelAerCheck = changeRxnBounds(modelAerCheck,modelAerCheck.rxns(end),0.01,'b');
    optSink = optimizeCbModel(modelAerCheck);
    growthCheck{i,1} = metsFound{i,1};
    growthCheck{i,2} = optSink.f;
end

%all of them can be produced, so lets add the new mets
phbIdx = findMetIDs(iDT1278,metsNotFound(:,1));
metsNotFound(:,2) = iDT1278.metNames(phbIdx);
metsNotFound(:,3) = iDT1278.metFormulas(phbIdx);
metsNotFound(:,4) = {0;0};

rxnsAdd = findRxnsFromMets(iDT1278,metsNotFound(:,1));
rxnsAdd(:,2) = iDT1278.rxnNames(findRxnIDs(iDT1278,rxnsAdd(:,1)));
rxnsAdd(:,3) = printRxnFormula(iDT1278,rxnsAdd(:,1));

modelPHB = addMultipleMetabolites(modelAer2,metsNotFound(:,1),'metNames',metsNotFound(:,2),'metFormulas',metsNotFound(:,3),'metCharges',metsNotFound(:,4));
for i = 1: size(rxnsAdd,1)
   modelPHB = addReaction(modelPHB,rxnsAdd{i,1},'reactionName',rxnsAdd{i,2},'reactionFormula',rxnsAdd{i,3});
   growthVal = optimizeCbModel(modelPHB);
   growthVal.f
end

[grow1,grow2] = exchangeSingleModel(modelPHB);



%compare with Thauera_aminoaromatica_S2.xml

thaueraS2 = readSBML('Thauera_aminoaromatica_S2.xml',1000);
rxnsRed = cellfun(@(rxns) rxns(3:end),thaueraS2.rxns,'UniformOutput',false);
thaueraS2.rxns = rxnsRed;
[thau1,thau2] = exchangeSingleModel(thaueraS2);
consComp = find(thau2.Flux<0);
redThau = thau2(consComp,:);
%exopolysacharides
%dtdprmn_c

exoPol = {'dtdprmn_c';'dtdp4d6dm_c';'dtdp4d6dg_c';'dtdpglu_c';'dtdp4addg_c';...
    'uacgam_c';'udpacgal_c';'uacmam_c';'uaccg_c';'unaga_c';'uaagmda_c';...
    'u3hga_c';'u3aga_c';'u23ga_c'};


%check first if we can use acetate

modelAc = changeRxnBounds(modelPHB,'EX_glc__D_e',0,'l');
modelAc = changeRxnBounds(modelAc,'EX_ac_e',-3,'l');
modelAc = changeRxnBounds(modelAc,'EX_o2_e',0,'l');
optimizeCbModel(modelAc)
[ac1,ac2] = exchangeSingleModel(modelAc);


%lets check if we can produce all the EPS and PHB with this set up

for i = 1 : length(exoPol)
    modelTest2 = addSinkReactions(modelAc,exoPol{i,1});
    modelTest2 = changeRxnBounds(modelTest2,modelTest2.rxns(end),0.01,'b');
    optSink = optimizeCbModel(modelTest2);
    growthCheck2{i,1} = exoPol{i,1};
    growthCheck2{i,2} = optSink.f;

end

%all of them can be produced in the presence of acetate!!

phbRules = {'PHBS_syn'};
phbRules(:,2) = {'(TMZ1T_RS00880) or (TMZ1T_RS05010) or (TMZ1T_RS05610)'};

findPHB = findRxnIDs(modelAc,phbRules(:,1));
modelAc.grRules{findPHB,1} = phbRules{1,2};
modelAc = updateGenes(modelAc);
modelAc = removeUnusedGenes(modelAc);
%add 34dhbz_e
metsArom = {'34dhbz_e';'4cml_c';'5odhf2a_c'};
metsArom(:,2) = {'3,4-Dihydroxybenzoate';'4-Carboxymuconolactone';'5-Oxo-4,5-dihydrofuran-2-acetate'};
metsArom(:,3) = {'C7H5O4';'C7H4O6';'C6H5O4'};

rxnsArom = {'3_4DHBZt2';'MUCCY_kt';'4CMLCL_kt';'OXOAEL'};
rxnsArom(:,2) = {'Protocatechuate outher membrane transport';'3 carboxy cis cis muconate cycloisomerase';...
    '4-CARBOXYMUCONOLACTONE-DECARBOXYLASE-RXN';'3-oxoadipate enol-lactone hydrolase'};
rxnsArom(:,3) = {'34dhbz_e <=> 34dhbz_p';'h_c + CCbuttc_c => 4cml_c';...
    'h_c + 4cml_c => co2_c + 5odhf2a_c';'h2o_c + 5odhf2a_c => h_c + 3oxoadp_c'};
rxnsArom(:,4) = {'';'TMZ1T_RS07310';'TMZ1T_RS07310';''};

modelArom = addMultipleMetabolites(modelAc,metsArom(:,1),'metNames',metsArom(:,2),'metFormulas',metsArom(:,3));
modelArom = addExchangeRxn(modelArom,{'34dhbz_e'},0,1000);
for i = 1: size(rxnsArom,1)
   modelArom = addReaction(modelArom,rxnsArom{i,1},'reactionName',rxnsArom{i,2},'reactionFormula',rxnsArom{i,3},'geneRule',rxnsArom{i,4});
   growthVal = optimizeCbModel(modelArom);
   growthVal.f
end

tableComp = readtable('carbonAndNitrogenThauera2.xlsx','Sheet',1);
compID = tableComp.IDCompound;
idxComp = find(~cellfun(@isempty,compID));
realComp = compID(idxComp);
realCompc = strcat(erase(realComp,'_e'),'_c');
realCompp = strcat(erase(realComp,'_e'),'_p');
checkArome = findMetIDs(modelArom,realComp);
checkAromc = findMetIDs(modelArom,realCompc);
checkAromp = findMetIDs(modelArom,realCompp);

%cytosol

cytFound = realCompc(checkAromc~=0);
cytNotFound = realCompc(checkAromc==0);

%exchange 
exFound = realComp(checkArome~=0);
exNotFound = realComp(checkArome==0);

%periplasm

perFound = realCompp(checkAromp~=0);
perNotFound = realCompp(checkAromp==0);


%new reactions to remove
aromRmv = {'4MPDH';'BZt'};

modelArom2 = removeRxns(modelArom,aromRmv);
compData1 = readtable('carbonAndNitrogenThauera4.xlsx','Sheet',1);
compData2 = readtable('carbonAndNitrogenThauera4.xlsx','Sheet',2);
compData3 = readtable('carbonAndNitrogenThauera4.xlsx','Sheet',3);

modelArom2 = addMultipleMetabolites(modelArom2,compData2.compound,'metNames',compData2.name,'metFormulas',compData2.formula);
for i = 1: size(compData3.reactionID,1)
   modelArom2 = addReaction(modelArom2,compData3.reactionID{i,1},'reactionName',compData3.reactionName{i,1},'reactionFormula',compData3.formula{i,1},'geneRule',compData3.genes{i,1});
   growthVal = optimizeCbModel(modelArom2);
   growthVal.f
end

modelArom2 = removeUnusedGenes(modelArom2);

aromExch = compData1.IDCompound;
aromRequired = compData1.requiredInModel;
exchAromRed = aromExch(aromRequired==1);
condExp = compData1.Condition(aromRequired==1);

[arom1,arom2] = exchangeSingleModel(modelArom2);
modExchArom = strcat('EX_',exchAromRed);
modelArom2 = changeRxnBounds(modelArom2,'BZt1pp',-1000,'l');

for i = 1 : length (modExchArom)
    condStr = 'Denitrification';
    compStr = strcmp(condExp{i,1},condStr);
    if compStr == 1 
        modelDent = changeRxnBounds(modelArom2,'EX_ac_e',0,'l');
        modelDent = changeRxnBounds(modelDent,'EX_o2_e',-2,'l');
        modelDent = changeRxnBounds(modelDent,modExchArom{i,1},-3,'l');
        modelDent = changeRxnBounds(modelDent,'EX_no3_e',-5,'l');
        optCond = optimizeCbModel(modelDent);
        %[check1,check2] = exchangeSingleModel(modelDent);
        expCheck{i,1} = optCond.f; 
        expCheck{i,2} = condExp{i,1};
        expCheck{i,3} = modExchArom{i,1};
    else
        modelArg = changeRxnBounds(modelArom2,'EX_ac_e',0,'l');
        modelArg = changeRxnBounds(modelArg,'EX_o2_e',-3,'l');
        modelArg = changeRxnBounds(modelArg,'EX_no3_e',0,'l');
        modelArg = changeRxnBounds(modelArg,modExchArom{i,1},-3,'l');
        optCond2 = optimizeCbModel(modelArg);
        %[check1,check2] = exchangeSingleModel(modelDent);
        expCheck{i,1} = optCond2.f;
        expCheck{i,2} = condExp{i,1};
        expCheck{i,3} = modExchArom{i,1};

    end
end

model4hbz = changeRxnBounds(modelArom2,'EX_ac_e',0,'l');
%fix benzene
checkBze = findRxnsFromMets(modelDent,'bz_c');
bzeForm = printRxnFormula(modelDent,checkBze);

%check
checkTest = findRxnsFromMets(modelDent,'for_c');
testForm = printRxnFormula(modelDent,checkTest);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%BIOLOG PLATES%%%%%%%%%%%%
%bring PM2Thauera sheet

PMData = readtable('biologInfo.xlsx','Sheet','PM2Thauera');

%aerobically grown in ammonium 
PMSubs = PMData.Substrates; PMbigg = PMData.ID; PMusage = PMData.Usage;
exIdx = find(contains(PMbigg,'_e'));
redPM = PMSubs(exIdx); redPM(:,2) =PMbigg(exIdx);
redPM(:,3) = PMusage(exIdx);
for i = 1 : length(redPM)
    checkExch = findRxnIDs(modelArom2,strcat('EX_',redPM{i,2}))==0;
    PMInfo{i,1} = redPM{i,2};
    if checkExch
        PMInfo{i,2} = NaN;
    else
        modelPM = changeRxnBounds(modelArom2,'EX_o2_e',-5,'l');
        modelPM = changeRxnBounds(modelPM,'EX_ac_e',0,'l');
        modelPM = changeRxnBounds(modelPM,'EX_no3_e',0,'l');
        modelPM = changeRxnBounds(modelPM,strcat('EX_',redPM{i,2}),-5,'l');
        growthPM = optimizeCbModel(modelPM);
        PMInfo{i,2} = growthPM.f;
    end
end

PMnum = cell2mat(PMInfo(:,2));
GrowthExp = PMInfo(PMnum>0);
nonGrowthExp = setdiff(PMInfo(:,1),GrowthExp);

%%%%%
%check which ones are in the model
%nonGrowthExp
exchIden = findRxnIDs(modelArom2,strcat('EX_',nonGrowthExp));
presentExch = nonGrowthExp(exchIden~=0);
missExch = nonGrowthExp(exchIden==0);

%%%
%lets work with presentExch
%fixExch = {};



modelFixExch = modelArom2;

for i = 1 : size(fixExch,1)
   modelFixExch = addReaction(modelFixExch,fixExch{i,1},'reactionName',fixExch{i,2},'reactionFormula',fixExch{i,3});
   growthVal = optimizeCbModel(modelFixExch);
   growthVal.f
end

[checkExch1, checkExch2] = exchangeSingleModel(modelFixExch);
%etha_e

modelEtha = changeRxnBounds(modelFixExch,'EX_etha_e',-5,'l');
modelEtha = changeRxnBounds(modelEtha,'EX_o2_e',-3,'l');
modelEtha = changeRxnBounds(modelEtha,'EX_ac_e',0,'l');
modelEtha = changeRxnBounds(modelEtha,'EX_no3_e',0,'l');
optimizeCbModel(modelEtha) %FIXED!


%glcn_e
%fixExch(2,:) = actualAdded(1,:);
modelglcn = changeRxnBounds(modelFixExch,'EX_glcn_e',-5,'l');
modelglcn = changeRxnBounds(modelglcn,'EX_o2_e',-3,'l');
modelglcn = changeRxnBounds(modelglcn,'EX_ac_e',0,'l');
modelglcn = changeRxnBounds(modelglcn,'EX_no3_e',0,'l');
optimizeCbModel(modelglcn) %FIXED!

%his__L_e
modelhisLe = changeRxnBounds(modelFixExch,'EX_his__L_e',-5,'l');
modelhisLe = changeRxnBounds(modelhisLe,'EX_o2_e',-3,'l');
modelhisLe = changeRxnBounds(modelhisLe,'EX_ac_e',0,'l');
modelhisLe = changeRxnBounds(modelhisLe,'EX_no3_e',0,'l');
optimizeCbModel(modelhisLe) %FIXED!

%peamn_e
modelpeamne = changeRxnBounds(modelFixExch,'EX_peamn_e',-5,'l');
modelpeamne = changeRxnBounds(modelpeamne,'EX_o2_e',-3,'l');
modelpeamne = changeRxnBounds(modelpeamne,'EX_ac_e',0,'l');
modelpeamne = changeRxnBounds(modelpeamne,'EX_no3_e',0,'l');
optimizeCbModel(modelpeamne) %FIXED!

%phe__L_e
modelpheLe = changeRxnBounds(modelFixExch,'EX_phe__L_e',-5,'l');
modelpheLe = changeRxnBounds(modelpheLe,'EX_o2_e',-3,'l');
modelpheLe = changeRxnBounds(modelpheLe,'EX_ac_e',0,'l');
modelpheLe = changeRxnBounds(modelpheLe,'EX_no3_e',0,'l');
optimizeCbModel(modelpheLe) %FIXED!

%thymd_e
modelthymde = changeRxnBounds(modelFixExch,'EX_thymd_e',-5,'l');
modelthymde = changeRxnBounds(modelthymde,'EX_o2_e',-3,'l');
modelthymde = changeRxnBounds(modelthymde,'EX_ac_e',0,'l');
modelthymde = changeRxnBounds(modelthymde,'EX_no3_e',0,'l');
optimizeCbModel(modelthymde) %FIXED!

%extraMets = {};
%extraRxns = {};

modelFixExch2=modelFixExch;
exchData1 = readtable('exchangeExtra.xlsx','Sheet',1);
exchData2 = readtable('exchangeExtra.xlsx','Sheet',2);
modelFixExch2 = addMultipleMetabolites(modelFixExch2,exchData2.mets,'metNames',exchData2.metNames,'metFormulas',exchData2.formulas);

for i = 1 : size(exchData1.rxns,1)
   modelFixExch2 = addReaction(modelFixExch2,exchData1.rxns{i,1},'reactionName',exchData1.rxnNames{i,1},'reactionFormula',exchData1.formula{i,1});
   growthVal = optimizeCbModel(modelFixExch2);
   growthVal.f
end
[exch1, exch2] = exchangeSingleModel(modelFixExch2);

%acon_C_e
modelacone = changeRxnBounds(modelFixExch2,'EX_acon_C_e',-5,'l');
modelacone = changeRxnBounds(modelacone,'EX_o2_e',-3,'l');
modelacone = changeRxnBounds(modelacone,'EX_ac_e',0,'l');
modelacone = changeRxnBounds(modelacone,'EX_no3_e',0,'l');
optimizeCbModel(modelacone) %FIXED!
[exch1, exch2] = exchangeSingleModel(modelacone);

%balamd_e
modelbalamde = changeRxnBounds(modelFixExch2,'EX_balamd_e',-5,'l');
modelbalamde = changeRxnBounds(modelbalamde,'EX_o2_e',-3,'l');
modelbalamde = changeRxnBounds(modelbalamde,'EX_ac_e',0,'l');
modelbalamde = changeRxnBounds(modelbalamde,'EX_no3_e',0,'l');
optimizeCbModel(modelbalamde) %FIXED!
[exch1, exch2] = exchangeSingleModel(modelbalamde);

%bhb_e
modelbhbe = changeRxnBounds(modelFixExch2,'EX_bhb_e',-5,'l');
modelbhbe = changeRxnBounds(modelbhbe,'EX_o2_e',-3,'l');
modelbhbe = changeRxnBounds(modelbhbe,'EX_ac_e',0,'l');
modelbhbe = changeRxnBounds(modelbhbe,'EX_no3_e',0,'l');
optimizeCbModel(modelbhbe) %FIXED!
[exch1, exch2] = exchangeSingleModel(modelbhbe);

%guln__L_e
modelgulnLe = changeRxnBounds(modelFixExch2,'EX_guln__L_e',-5,'l');
modelgulnLe = changeRxnBounds(modelgulnLe,'EX_o2_e',-3,'l');
modelgulnLe = changeRxnBounds(modelgulnLe,'EX_ac_e',0,'l');
modelgulnLe = changeRxnBounds(modelgulnLe,'EX_no3_e',0,'l');
optimizeCbModel(modelgulnLe) %FIXED!
[exch1, exch2] = exchangeSingleModel(modelgulnLe);

%urcan_e
modelurcane = changeRxnBounds(modelFixExch2,'EX_urcan_e',-5,'l');
modelurcane = changeRxnBounds(modelurcane,'EX_o2_e',-3,'l');
modelurcane = changeRxnBounds(modelurcane,'EX_ac_e',0,'l');
modelurcane = changeRxnBounds(modelurcane,'EX_no3_e',0,'l');
optimizeCbModel(modelurcane) %FIXED!

%btd_RR_e


%lets run again the PM analysis 
%lets determine the usage of the compounds
redPMMod = redPM;
redPMMod{17,2} = 'gln__L_e';
redPMMod(16,:) = [];
clear PMInfo2
for i = 1 : length(redPMMod)
    checkExch = findRxnIDs(modelFixExch2,strcat('EX_',redPMMod{i,2}))==0;
    PMInfo2{i,1} = redPMMod{i,2};
    if checkExch
        PMInfo2{i,2} = NaN;
    else
        modelPM = changeRxnBounds(modelFixExch2,'EX_o2_e',-5,'l');
        modelPM = changeRxnBounds(modelPM,'EX_ac_e',0,'l');
        modelPM = changeRxnBounds(modelPM,'EX_no3_e',0,'l');
        modelPM = changeRxnBounds(modelPM,strcat('EX_',redPMMod{i,2}),-4,'l');
        growthPM = optimizeCbModel(modelPM);
        PMInfo2{i,2} = growthPM.f;
    end
end

%PM final data

dataPM = redPMMod;
dataPM(:,4) = PMInfo2(:,2);
growVal = cell2mat(dataPM(:,4));
growVal(22,1) = 0;
growVal(29,1) = 0;

%lets divide in high, medium, and low growth

highG = contains(dataPM(:,3),'high');
mediumG = contains(dataPM(:,3),'medium');
lowG = contains(dataPM(:,3),'low');
highData = dataPM(highG,:);mediumData = dataPM(mediumG,:);lowData = dataPM(lowG,:);
avgHigh = mean(growVal(highG)); avgMed = mean(growVal(mediumG));
avgLow = mean(growVal(lowG));
stdHigh = std(growVal(highG)); stdMed = std(growVal(mediumG));
stdLow = std(growVal(lowG));
[sortData, sortIdx] = sort(growVal,'descend');
reorgData = dataPM(sortIdx,:);

modelFixExch3 = updateGenes(modelFixExch2);
modelFixExch3 = removeUnusedGenes(modelFixExch3);
modelFixExch3 = generateRules(modelFixExch3);
modelFixExch3 = updateGenes(modelFixExch3);
modelFixExch3 = removeUnusedGenes(modelFixExch3);
[exchD1,exchD2] = exchangeSingleModel(modelFixExch3);

%

%lets annotate all the subsystems for every reaction
thaueraSubs1 = readcell('modelSubsThauera.xlsx','Sheet',2);
rulesSubs = thaueraSubs1(:,2);
for i = 1:numel(rulesSubs)
    if ismissing(rulesSubs{i})
        rulesSubs{i} = [];
    end
end
rulesSubs = cellfun(@char,rulesSubs,'UniformOutput',false);
emptySubs = find(cellfun(@isempty,rulesSubs));
rxnsSubs = thaueraSubs1(:,1);
rxnsEmpty = rxnsSubs(emptySubs);
checkRxns2 = findRxnIDs(iML1515,rxnsEmpty);
redID = checkRxns2(checkRxns2~=0);
rxnsIden = rxnsEmpty(checkRxns2~=0);
posRxns = cell2mat(cellfun(@(pos) strmatch(pos,rxnsSubs,'exact'),rxnsIden,'UniformOutput',false));
thaueraSubs1(posRxns,2) = iML1515.subSystems(redID);

for i = 1:length(thaueraSubs1)
    if ismissing(thaueraSubs1{i,2})
        thaueraSubs1{i,2} = [];
    end
end

thaueraSubs1(:,1) = cellfun(@char,thaueraSubs1(:,1),'UniformOutput',false);
thaueraSubs1(:,2) = cellfun(@char,thaueraSubs1(:,2),'UniformOutput',false);



%single gene deletion
clear grDataT, clear fluxGenesPM
for i = 1 : length(redPMMod)
    i
    clear grRatio, clear grRateKO, clear grRateWT, clear has Effect
    clear delRxns, clear fluxSolution
    checkGrowth = growVal(i)==0;
    if checkGrowth
        grDataT{i,1} = {}; grDataT{i,2} = {}; grDataT{i,3} = {};
        grDataT{i,4} = {}; grDataT{i,5} = {}; grDataT{i,6} = {};

    else
        modelPM = changeRxnBounds(modelFixExch3,'EX_o2_e',-5,'l');
        modelPM = changeRxnBounds(modelPM,'EX_ac_e',0,'l');
        modelPM = changeRxnBounds(modelPM,'EX_no3_e',0,'l');
        modelPM = changeRxnBounds(modelPM,strcat('EX_',redPMMod{i,2}),-4,'l');
        solGenes = optimizeCbModel(modelPM,'max')
        fluxGenesPM{i,1} = solGenes.x;
        [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution] = singleGeneDeletion2(modelPM);
        grRatio(isnan(grRatio)) = 0; grRateKO(isnan(grRateKO)) = 0;
        grDataT{i,1} = grRatio; grDataT{i,2} = grRateKO; grDataT{i,3} = grRateWT;
        grDataT{i,4} = hasEffect; grDataT{i,5} = delRxns; grDataT{i,6} = fluxSolution; 
    end
end

[check1,check2] = exchangeSingleModel(modelPM);
%lets remove the conditions without any growth
emptyData = find(~cellfun(@isempty,grDataT(:,1)));
fluxGenesRed = fluxGenesPM(emptyData);
fluxGenes = horzcat(fluxGenesRed{:});
grDataTR = grDataT(emptyData,:); ratioData = grDataTR(:,1);
ratioData = horzcat(ratioData{:});
ratio = cellfun(@sum, mat2cell(ratioData, ones(size(ratioData, 1), 1), size(ratioData, 2)));
nSize = size(grDataTR,1);
eqGrowth=find(ratio==nSize);affGrowth = find(ratio>0 & ratio<nSize);
affGrowthMod = find(ratio>0 & ratio<35);
zeroGrowth = find(ratio<=0);
genesNoImp = modelFixExch3.genes(eqGrowth);
genesSomeImp = modelFixExch3.genes(affGrowth);
genesSomeImp2 = modelFixExch3.genes(affGrowthMod);
lethalGenes = modelFixExch3.genes(zeroGrowth);
%genes analysis
%lethal genes
[~, ListResults] = findRxnsFromGenes(modelFixExch3,lethalGenes,'ListResultsFlag',1);
subsLethal = ListResults(:,3);
uniqLethal = unique(subsLethal);
lethalDistb = cell2mat(cellfun(@(str) length(strmatch(str,subsLethal,'exact')),uniqLethal,'UniformOutput',false));
lethalPerc = lethalDistb/sum(lethalDistb)*100;
pie(lethalPerc)
l = legend(uniqLethal);
hold on
title('Subsystems affected by lethal genes under aerobic conditions')
bussyMedium = redPMMod(emptyData,2); bussyMedium=erase(bussyMedium,'_e');

%genes which are lethal in only some conditions
ratioImp = ratioData(affGrowth,:);
%plot considering row standarization
cgo = clustergram(ratioImp,'Colormap',redbluecmap,'Standardize','Row')
set(cgo,'ColumnLabels',bussyMedium)
addTitle(cgo,'Genes affected for different carbon sources under aerobic conditions')
addXLabel(cgo,'Carbon source');addYLabel(cgo,'Genes')
plot(cgo)
%reduced genes which are lethal in only some conditions
ratioImp2 = ratioData(affGrowthMod,:);
cgo = clustergram(ratioImp2,'Colormap',redbluecmap,'Standardize','Row')
set(cgo,'RowLabels',erase(genesSomeImp2,'_'),'ColumnLabels',bussyMedium)
addTitle(cgo,'Genes affected for different carbon sources under aerobic conditions')
addXLabel(cgo,'Carbon source');addYLabel(cgo,'Genes')
plot(cgo)

%lets evaluate which growing conditions are capable to produce higher
%concentrations of PHBs
bussyMediumMod = strcat(bussyMedium,'_e');
phbIdx = findRxnIDs(modelFixExch3,'PHBS_syn');
clear maxProd
for i = 1 : length(bussyMedium)
    i
    modelPM = changeRxnBounds(modelFixExch3,'EX_o2_e',-5,'l');
    modelPM = changeRxnBounds(modelPM,'EX_ac_e',0,'l');
    modelPM = changeRxnBounds(modelPM,'EX_no3_e',0,'l');
    modelPM = changeRxnBounds(modelPM,strcat('EX_',bussyMediumMod{i,1}),-4,'l');
    growthPM = optimizeCbModel(modelPM);
    flux = 0;
    while growthPM.f>0 || ~isnan(growthPM.f)
        flux = flux + 0.001;
        modelPM2 = changeRxnBounds(modelPM,'PHBS_syn',flux,'b');
        growthPM = optimizeCbModel(modelPM2);
    end
    modelPM2 = changeRxnBounds(modelPM,'PHBS_syn',flux-0.001,'b');
    growthPM = optimizeCbModel(modelPM2);
    maxProd(i,1) = growthPM.f;  maxProd(i,2) = growthPM.x(phbIdx);
    
end
%exopolymer precursors exoPol

rxnList = cellfun(@(mets) findRxnsFromMets(modelFixExch3, mets),exoPol,'UniformOutput', false);
clear growthProd2, clear maxProd2
for i = 1 : length(bussyMedium)
    i
    modelPM = changeRxnBounds(modelFixExch3,'EX_o2_e',-5,'l');
    modelPM = changeRxnBounds(modelPM,'EX_ac_e',0,'l');
    modelPM = changeRxnBounds(modelPM,'EX_no3_e',0,'l');
    modelPM = changeRxnBounds(modelPM,strcat('EX_',bussyMediumMod{i,1}),-4,'l');
    growthPM = optimizeCbModel(modelPM);
    for j = 1 : length(exoPol)
        flux = 0;
        modelPM2 = addSinkReactions(modelPM,exoPol{j,1});
        modelPM2 = changeRxnBounds(modelPM2,modelPM2.rxns(end),flux,'b');
        growthPM = optimizeCbModel(modelPM2);
        if growthPM.f>0 || ~isnan(growthPM.f)
            while growthPM.f>0 || ~isnan(growthPM.f)
                flux = flux + 0.005; 
                modelPM2 = changeRxnBounds(modelPM2,modelPM2.rxns(end),flux,'b');
                growthPM = optimizeCbModel(modelPM2);

            end
            modelPM2 = changeRxnBounds(modelPM2,modelPM2.rxns(end),flux-0.005,'b');
            growthPM = optimizeCbModel(modelPM2);
            growthProd2(i,j) = growthPM.f;
            maxProd2(i,j) = growthPM.x(end);
        else
            growthProd2(i,j) = 0;
            maxProd2(i,j) = 0;
        end
        
    end
    
end


%exopolymer precursors exoPol using concentration g/L
MW = computeMW(modelFixExch3,bussyMediumMod);

rxnList = cellfun(@(mets) findRxnsFromMets(modelFixExch3, mets),exoPol,'UniformOutput', false);
clear growthProd4, clear maxProd4
for i = 1 : length(bussyMedium)
    i
    modelPM = changeRxnBounds(modelFixExch3,'EX_o2_e',-5,'l');
    modelPM = changeRxnBounds(modelPM,'EX_ac_e',0,'l');
    modelPM = changeRxnBounds(modelPM,'EX_no3_e',0,'l');
    modelPM = changeRxnBounds(modelPM,strcat('EX_',bussyMediumMod{i,1}),-100/MW(i),'l');
    growthPM = optimizeCbModel(modelPM);
    for j = 1 : length(exoPol)
        flux = 0;
        modelPM2 = addSinkReactions(modelPM,exoPol{j,1});
        modelPM2 = changeRxnBounds(modelPM2,modelPM2.rxns(end),flux,'b');
        growthPM = optimizeCbModel(modelPM2);
        if growthPM.f>0 || ~isnan(growthPM.f)
            while growthPM.f>0 || ~isnan(growthPM.f)
                flux = flux + 0.003; 
                modelPM2 = changeRxnBounds(modelPM2,modelPM2.rxns(end),flux,'b');
                growthPM = optimizeCbModel(modelPM2);

            end
            modelPM2 = changeRxnBounds(modelPM2,modelPM2.rxns(end),flux-0.003,'b');
            growthPM = optimizeCbModel(modelPM2);
            growthProd4(i,j) = growthPM.f;
            maxProd4(i,j) = growthPM.x(end);
        else
            growthProd4(i,j) = 0;
            maxProd4(i,j) = 0;
        end
        
    end
    
end

%Plotting exopolymers mass concentration
cgo4 = clustergram(maxProd4','Colormap',redbluecmap,'Standardize','Row')
set(cgo4,'RowLabels',exoPol,'ColumnLabels',bussyMedium)
addTitle(cgo4,'Exopolmyer production fluxes for different carbon sources under aerobic conditions')
addXLabel(cgo4,'Carbon source');addYLabel(cgo4,'Exopolymers')
plot(cgo4)

cgo5 = clustergram(maxProd4','Colormap',redbluecmap,'Standardize','Column')
set(cgo5,'RowLabels',exoPol,'ColumnLabels',bussyMedium)
addTitle(cgo5,'Exopolmyer production fluxes for different carbon sources under aerobic conditions')
addXLabel(cgo5,'Carbon source');addYLabel(cgo5,'Exopolymers')


%Bring subsystems curated information


thaueraSubs2 = readcell('modelSubsThauera.xlsx','Sheet',2);
modelFixExch3.subSystems = thaueraSubs2(:,2);
uniqSubs = unique(modelFixExch3.subSystems);
accumSubs = zeros(length(uniqSubs),size(fluxGenes,2));
for i = 1 : size(fluxGenes,2)
    clear sumSubs
    actIdx = find(fluxGenes(:,i)~=0);
    actFluxes = fluxGenes(actIdx,i);
    actRxns = modelFixExch3.rxns(actIdx);
    actSubs = modelFixExch3.subSystems(actIdx);
    uniqAct = unique(actSubs);
    subsIdx = cell2mat(cellfun(@(str) strmatch(str,uniqSubs,'exact'),uniqAct,'UniformOutput',false));
    posSubs = cellfun(@(str) strmatch(str,actSubs,'exact'),uniqAct,'UniformOutput',false);
    for j = 1 : length(posSubs)
        sumSubs(j,1) = sum(abs(actFluxes(posSubs{j,1})));  
    end
    accumSubs(subsIdx,i) = sumSubs;
end
ratioSum = cellfun(@sum, mat2cell(accumSubs, ones(size(accumSubs, 1), 1), size(accumSubs, 2)));
bussySubs = find(ratioSum~=0);
redAccumSubs = accumSubs(bussySubs,:);
redSubs = uniqSubs(bussySubs);
notUsedSubs = uniqSubs(ratioSum==0);
%plot considering row standarization
cgo = clustergram(redAccumSubs,'Colormap',redbluecmap,'Standardize','Row')
set(cgo,'RowLabels',redSubs,'ColumnLabels',bussyMedium)
addTitle(cgo,'Subsystem flux distributions for different carbon sources under aerobic conditions')
addXLabel(cgo,'Carbon source');addYLabel(cgo,'Subsystem')

%plot considering column standarization
cgo2 = clustergram(redAccumSubs,'Colormap',redbluecmap,'Standardize','Column')
set(cgo2,'RowLabels',redSubs,'ColumnLabels',bussyMedium)
addTitle(cgo2,'Subsystem flux distributions for different carbon sources under aerobic conditions')
addXLabel(cgo2,'Carbon source');addYLabel(cgo2,'Subsystem')

%single reaction deletion
clear grDataT2
for i = 1 : length(redPMMod)
    i
    clear grRatio, clear grRateKO, clear grRateWT, clear has Effect
    clear delRxns, clear fluxSolution
    checkGrowth = growVal(i)==0;
    if checkGrowth
        grDataT2{i,1} = {}; grDataT2{i,2} = {}; grDataT2{i,3} = {};
        grDataT2{i,4} = {};

    else
        modelPM = changeRxnBounds(modelFixExch3,'EX_o2_e',-5,'l');
        modelPM = changeRxnBounds(modelPM,'EX_ac_e',0,'l');
        modelPM = changeRxnBounds(modelPM,'EX_no3_e',0,'l');
        modelPM = changeRxnBounds(modelPM,strcat('EX_',redPM{i,2}),-4,'l');
        grRateWT=optimizeCbModel(modelPM,'max');
        showprogress(0,'Single gene deletion analysis in progress ...');
        for j = 1 : length(modelPM.rxns)
            showprogress(j/length(modelPM.rxns));
            modelPM2 = changeRxnBounds(modelPM,modelPM.rxns{j,1},0,'b');
            grRateKO = optimizeCbModel(modelPM2,'max');
            if isnan(grRateKO.f)
                grRateKO2(j,1) = 0; 
                grRatio(j,1) = 0 ;
                fluxSolution(:,j) = zeros(length(modelPM.rxns),1);
            else
                grRateKO2(j,1) = grRateKO.f; 
                grRatio(j,1) = grRateKO.f/grRateWT.f ;
                fluxSolution(:,j) = grRateKO.x;
            end
            
            
            
        end
        grRatio(isnan(grRatio)) = 0;
        grDataT2{i,1} = grRatio; grDataT2{i,2} = grRateKO2; grDataT2{i,3} = grRateWT.f;
        grDataT2{i,4} = fluxSolution; 
    end
end

%finding lethal reactions, affect reactions and non-relevant reactions

emptyData2 = find(~cellfun(@isempty,grDataT2(:,1)));
grDataTR2 = grDataT2(emptyData2,:); ratioData2 = grDataTR2(:,1);
ratioData2 = horzcat(ratioData2{:});
ratio2 = cellfun(@sum, mat2cell(ratioData2, ones(size(ratioData2, 1), 1), size(ratioData2, 2)));
nSize2 = size(grDataTR2,1);
eqGrowth2=find(ratio2==nSize2);affGrowth2 = find(ratio2>0 & ratio2<nSize2);
zeroGrowth2 = find(ratio2<=0);
rxnsNoImp = modelFixExch3.rxns(eqGrowth2);
rxnsSomeImp = modelFixExch3.rxns(affGrowth2);
lethalRxns = modelFixExch3.rxns(zeroGrowth2);

%check EPS production
for i = 1 : length(exoPol)
    modelTest2 = addSinkReactions(modelPM,exoPol{i,1});
    modelTest2 = changeRxnBounds(modelTest2,modelTest2.rxns(end),0.01,'b');
    optSink = optimizeCbModel(modelTest2);
    growthCheck3{i,1} = exoPol{i,1};
    growthCheck3{i,2} = optSink.f;

end

[exchC,exchD] = exchangeSingleModel(modelFixExch3);
uniqSubs3 = unique(modelFixExch3.subSystems);
subsIdxT = cell2mat(cellfun(@(str) length(strmatch(str,modelFixExch3.subSystems,'exact')),uniqSubs3,'UniformOutput',false));

%share iYL1228 and iDT863
%share1228 = {};

%lets try setting one concentration 10g/L
%lets evaluate which growing conditions are capable to produce higher
%concentrations of PHBs
bussyMediumMod = strcat(bussyMedium,'_e');
phbIdx = findRxnIDs(modelFixExch3,'PHBS_syn');
clear maxProd3
for i = 1 : length(bussyMedium)
    i
    modelPM = changeRxnBounds(modelFixExch3,'EX_o2_e',-5,'l');
    modelPM = changeRxnBounds(modelPM,'EX_ac_e',0,'l');
    modelPM = changeRxnBounds(modelPM,'EX_no3_e',0,'l');
    modelPM = changeRxnBounds(modelPM,strcat('EX_',bussyMediumMod{i,1}),-100/MW(i),'l');
    growthPM = optimizeCbModel(modelPM);
    flux = 0;
    while growthPM.f>0 || ~isnan(growthPM.f)
        flux = flux + 0.001;
        modelPM2 = changeRxnBounds(modelPM,'PHBS_syn',flux,'b');
        growthPM = optimizeCbModel(modelPM2);
    end
    modelPM2 = changeRxnBounds(modelPM,'PHBS_syn',flux-0.001,'b');
    growthPM = optimizeCbModel(modelPM2);
    maxProd3(i,1) = growthPM.f;  maxProd3(i,2) = growthPM.x(phbIdx);
    
end
maxProd4 = maxProd3;
maxProd4(:,2) = maxProd4(:,2)./max(maxProd4(:,2));
[~,sortIdx] = sort(maxProd4(:,2));
listPHB = horzcat(bussyMedium(sortIdx),num2cell(maxProd4(sortIdx,2)));

epsNames = modelFixExch3.metNames(findMetIDs(modelFixExch3, exoPol));
iDT863 = modelFixExch3;