%Made by Paweł Kruczkiewicz, 09-01-2021%

function u = skryptMES(n)
    %Gaussian 2 point quadrature for aproximating integral value% 
    function res = gaussian2pointQuad(f, a, b)
        c = (b-a)/2;
        d = (a+b)/2;
        epsilon = [-1/sqrt(3), 1/sqrt(3)];
        
        x = c.*epsilon + d;
        y = f(x);
        res = c*sum(y);
    end

    function [a, b] = computeIntervalToIntegrate(i, j)
        if i == j
            a_ind = max(i-1,1);
            b_ind = min(i+1,n+1);
            a = t(a_ind);
            b = t(b_ind);
        elseif abs(i - j) == 1
            if i > j
                [j, i] = deal(i, j);
            end
            a = t(i);
            b = t(j);
        else
            a = 0;
            b = 0;
        end
    end

    %1-dimentional B-spline of i-th index 
    function res = e(i)
        component1 = @(x) 0;
        component2 = @(x) 0;
        if i ~= n+1
            component1 = @(x) ((x >= t(i)) & (x < t(i+1)))*(t(i+1) - x)/h;
        end
        if i ~= 1
            component2 = @(x) ((x >= t(i-1)) & (x < t(i)))*(x - t(i-1))/h;
        end
      
        res = @(x) component1(x) + component2(x);   
    end

    %1-dimentional B-spline's derivative of i-th index 
    function res = der_e(i)
        component1 = @(x) 0;
        component2 = @(x) 0;
        if i ~= n+1
            component1 = @(x) ((x >= t(i)) & (x < t(i+1)))*(-1)/h ;
        end
        if i ~= 1
            component2 = @(x) ((x >= t(i-1)) & (x < t(i)))*1/h;
        end
      
        res = @(x) component1(x) + component2(x);   
    end

    %value of function B (defined in derivation) on position (i,j)
    function res = B(i,j)
        ei = e(i);
        ej = e(j);
        der_ei = der_e(i);
        der_ej = der_e(j);
        
        f = @(x) k(x).*der_ei(x).*der_ej(x);
        [a,b] = computeIntervalToIntegrate(i, j);
        
        
        integral = gaussian2pointQuad(f, a, b);
        if a < 1 && 1 < b
            integral = gaussian2pointQuad(f, a, 1) + gaussian2pointQuad(f, 1, b);
        end
        
        res = integral - k(0).*ei(0).*ej(0);
    end

    function res = L(i)
        ei = e(i);
        res = -20.*k(0).*ei(0);
    end


    %k function - value 1 for x in [0,1], value 2 for x in (1,2]%
    k = @(x) (1*((x>=0) & (x <=1)) + 2*((x > 1) & (x <= 2)));
    
    %h - length of single interval%
    h = 2/n;
    
    %t - points of division
    t = linspace(0, 2, n+1);
    
    A = zeros(n+1);
    for i = 1:n
        A(i,i+1) = B(i+1,i);
        A(i,i) = B(i,i);
        A(i+1,i) = B(i,i+1); 
    end
    A(n+1, n) = 0;
    A(n, n+1) = 0;
    A(n+1, n+1) = 1;
    
    C = zeros(n+1, 1);
    for i = 1:(n+1)
        C(i, 1) = L(i);
    end
    
    u = A\C;
    
    plot(t,u, t, u, 'o');
end



