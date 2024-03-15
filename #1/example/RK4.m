%10阶4级强保稳，SSP-RK4格式
function Qn = RK4(A,Q,t)
Q1 = Q;
Q2 = Q;
for i = 1:5
Q1 = Q1 + (t/6)*L(A,Q1);
end
Q2 = (1/25)*Q2 + (9/25)*Q1;
Q1 = 15*Q2 - 5*Q1;
for i = 6:9
Q1 = Q1 + (t/6)*L(A,Q1);
end
Qn = Q2 + (3/5)*Q1 + (t/10)*L(A,Q1);
end