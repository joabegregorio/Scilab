//TARTAGLIA-CARDANO METHOD

function [raiz1,raiz2,raiz3] = Tartaglia_Cardano(a,b,c,d)

//turn the original equation into an x^3 + Ax^2+ Bx+ C
A = b / a;
B = c / a;
C = d / a;

//convertion constants to y^2 + py + q
p = B - (A^2)/3;
q = C + (2*(A^3))/27 - (A*B)/3;

//discriminant
delta =(q^2)/4 + (p^3)/27;

    if delta >= 0 then
        //first root
        y1 = nthroot((-q/2 + sqrt(delta)),3) + nthroot((-q/2 - sqrt(        delta)),3);
        raiz1 = y1- A/3;
        //secundary discriminant of an quadratic euqation
        delta2 = -3*(y1^2)- 4*p;
        //other 2 roots
        raiz2 = ((-y1 + sqrt(delta2))/2) -A/3;
        raiz3 = ((-y1 - sqrt(delta2))/2) - A/3;
    elseif delta <= 0 then 
        //find the roots by the Euler method
        ro   = sqrt(((q^2)/4)+ abs(delta));
        teta = acos(-q/(2*ro));
        raiz1 = (2*nthroot(ro,3)*cos(teta/3)) - A/3;
        raiz2 = (2*nthroot(ro,3)*cos((teta+(2*%pi))/3)) -A/3;
        raiz3 = (2*nthroot(ro,3)*cos((teta+(4*%pi))/3)) -A/3;
    end


endfunction
