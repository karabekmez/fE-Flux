%   fE-FLUX: feasible E-FLUX implementation
%   Created by Muhammed Erkan Karabekmez, PhD (22/09/2023)
%   Adjusts lb and ub of a model proportional to reaction expression values
%   within FVA limits.
%   Any adjustment that would failing growth simulation by FBA is not
%   allowed, in order to get feasible models.
%   Inputs: model, reaction_expression 
%   Output: bounded model
%   Dependencies: FBA and FVA functions (optimizeCbModel & fluxVariability in COBRA Toolbox)

model_RAW=model;                
Rxnexp=reaction_expression;     
n_Rxnexp=log2(Rxnexp);  % If transcriptomic data is already log2 normalized then by-pass this line       
[minFlux, maxFlux] = fluxVariability(model_RAW);        
for x=1:length(model_RAW.rxns)
    if (isnan(n_Rxnexp(x))==0 && Rxnexp(x)~=0) 
        model_RAW.lb(x)=((n_Rxnexp(x)/max(n_Rxnexp)))*minFlux(x);
        sol=optimizeCbModel(model_RAW);
        exist sol.f1;
        if (isequal(sol.origStat ,'INFEASIBLE') || (isnan(sol.f1)==1))
            model_RAW.lb(x)=minFlux(x);
        end
    end
    if (isnan(n_Rxnexp(x))==0 && Rxnexp(x)~=0)
        model_RAW.ub(x)=(n_Rxnexp(x)/max(n_Rxnexp))*maxFlux(x);
        sol=optimizeCbModel(model_RAW);
        exist sol.f1;
        if (isequal(sol.origStat ,'INFEASIBLE') || (isnan(sol.f1)==1))
            model_RAW.ub(x)=maxFlux(x);
        end
    end
end
save model_RAW model_RAW
