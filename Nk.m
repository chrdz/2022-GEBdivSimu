function res=Nk(x,k,Ne,he, xList)
% shape function Nk(x)
    if k == 1
%         f = interp1(ft,f,t); % Interpolate the data set (ft,f) at time t
%         g = interp1(gt,g,t); % Interpolate the data set (gt,g) at time t
%         dydt = -f.*y + g; % Evaluate ODE at time t
        indic = (x <= he);
        %res = (1-3*(x/he)+2*(x/he)^2)*indic;
        
        xi = x/(xList(3) - xList(1));
        res = Ntild1(xi)*indic;
        
    elseif k == 2*Ne+1
        %indic = (x >= (Ne-1)*he);
        %res = (-((x-(Ne-1)*he)/he)+2*((x-(Ne-1)*he)/he)^2)*indic;
        
        indic = (x >= xList(k-2));
        xi = (x - xList(k-2))/(xList(k) - xList(k-2)); 
        res = Ntild3(xi)*indic;
    elseif mod(k,2) ~= 0
        %indic1 = (x >= (k-3)*he/2  &&  x <= (k-1)*he/2);
        %indic2 = (x > (k-1)*he/2  &&  x <= (k+1)*he/2);
        %res = (-((x-(k-3)*he/2)/he)+2*((x-(k-3)*he/2)/he)^2)*indic1 + ...
        %    (1-3*((x-(k-1)*he/2)/he)+2*((x-(k-1)*he/2)/he)^2)*indic2;
        
        indic1 = (x >= xList(k-2)  &&  x <= xList(k));
        indic2 = (x > xList(k)  &&  x <= xList(k+2));
        xi1 = (x - xList(k-2))/(xList(k) - xList(k-2));
        xi2 = (x - xList(k))/(xList(k+2) - xList(k));
        res = Ntild3(xi1)*indic1 + Ntild1(xi2)*indic2;
        
    else
        %indic = (x >= (k-2)*he/2  &&  x <= k*he/2);
        %res = 4*(((x-(k-2)*he/2)/he)-((x-(k-2)*he/2)/he)^2)*indic; 
        
        indic = (x >= xList(k-1)  &&  x <= xList(k+1));
        xi = (x-xList(k-1))/(xList(k+1)-xList(k-1));
        res = Ntild2(xi)*indic;
    end
end