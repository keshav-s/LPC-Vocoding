function [a, k] = levinson792(phi)
% Performs Levinson-Derbin recursion
P = length(phi)-1;
E = phi(1);

k = zeros(P,1);
a = zeros(P,1);

k(1) = phi(2)/E;
a(1) = k(1);
prev_a = a;
E = (1-k(1)^2) * E;

for i=2:P
    idxs = 1:i-1;
    idxs = i+1 - idxs;
    k(i) = ( phi(i+1) - sum(prev_a(1:i-1).* phi(idxs)) ) / E;

    a(i) = k(i);
    for j=1:i-1
        a(j) = prev_a(j) - k(i)*prev_a(i-j);
    end

    E = (1-k(i)^2) * E;
    prev_a = a;
end

end
