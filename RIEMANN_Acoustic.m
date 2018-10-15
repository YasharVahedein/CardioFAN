function [PP,UU]=RIEMANN_Acoustic(P1,P2,U1,U2,roc)
PP=0.5*(P1+P2)-0.5*roc*(U2-U1);
UU=0.5*(U1+U2)-0.5*(P2-P1)/roc;
end

