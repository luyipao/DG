[X,T,Q,C] = LegendRKDG(0,2*pi,4*pi,3,4,100,0.1);
M = size(Q,2);
figure
xlabel('X');
ylabel('Q');
s = [X(1),X(end),-2,2];
axis(s);
[p,q] = size(Q);
Q2 = zeros(p,q);
for i = 1:p
    for j = 1:q
        Q2(i,j) = sin(X(i) - T(j));
    end
end
for i = 1:10:M
   plot(X,Q(:,i),X,Q2(:,i),'r.');
   axis(s);
   pause(0.001);
end