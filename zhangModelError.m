function [objectiveVal] = zhangModelError(paramCalib, q, rf, pet, sf, nMonths, tmax, tmin)

[qEst evap stor snowFinal storFinal groundFinal Snowpack] = zhang_model_snow(rf,pet,sf,nMonths,paramCalib(:,1),paramCalib(:,2),...
    paramCalib(:,3),paramCalib(:,5),paramCalib(:,4),0,0,0, tmax, tmin);
[qEst evap stor snowFinal storFinal groundFinal Snowpack] = zhang_model_snow(rf,pet,sf,nMonths,paramCalib(:,1),paramCalib(:,2),...
    paramCalib(:,3),paramCalib(:,5),paramCalib(:,4),storFinal,groundFinal,snowFinal, tmax, tmin);

temp1 = 0;
temp2 = 0;
q_bar = mean(q);
for i2 = 1:length(qEst)
temp1 = temp1 +(log(qEst(i2)) - log(q(i2)))^2;
temp2 = temp2 +(log(q(i2))-log(q_bar))^2;
end
fun1 = temp1/temp2;

temp1 = 0;
temp2 = 0;
for i3 = 1:size(qEst)
temp1 = temp1 +(qEst(i3) - q(i3))^2;
temp2 = temp2 +(q(i3)- q_bar)^2;
end
fun2 = temp1/temp2;

temp1 = 0;
temp2 = 0;
temp3 = 0;
qEst_bar = mean(qEst);
for i4 = 1:size(qEst)
temp1 = temp1 +(qEst(i4) - qEst_bar)*(q(i4)-q_bar);
temp2 = temp2 +(qEst(i4)- qEst_bar)^2;
temp3 = temp3 + (q(i4) - q_bar)^2;
end
fun3 = temp1/sqrt(temp2*temp3);

fun4 = abs(log(sum(q)/sum(qEst)));

objectiveVal = mean([fun1 fun2 fun3 fun4]);
end