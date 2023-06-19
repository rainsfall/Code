ops = sdpsettings('solver', 'Gurobi+', 'verbose', 2, 'debug', 1, 'gurobi.NonConvex', 2);
obj = sum(P_itH(:,2));
result = optimize(Constraints,obj,ops);

[model,recoverymodel,diagnostic,internalmodel] = export(Constraints, obj,sdpsettings('solver','GUROBI+'));    
iis = gurobi_iis(model);
gurobi_write(model, 'TestModel.lp');
find(iis.Arows)-1