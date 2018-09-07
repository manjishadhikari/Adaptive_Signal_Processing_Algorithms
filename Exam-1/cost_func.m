function [J]=cost_func(R,s1,s2,w)
    J= w'*R*w+abs(w'*s1-1)^2+abs(w'*s2)^2;
return