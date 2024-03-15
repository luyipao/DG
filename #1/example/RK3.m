%3阶强保稳龙格库塔方法
function Qn = RK3(A,Q,t)
k1 = Q + t*L(A,Q);
k2 = (3/4)*Q + (1/4)*k1 + (1/4)*t*L(A,k1);
Qn = (1/3)*Q + (2/3)*k2 + (2/3)*t*L(A,k2);
end