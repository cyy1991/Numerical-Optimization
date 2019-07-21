% 
function [res_t, res_a] = temp(lam_his, m1)

    ReadData();
    
    g = @(lam) quadgk(@(x) exp(lam' * [x.^0;x.^1;x.^2;x.^3;x.^4;x.^5;x.^6;x.^7;x.^8;x.^9;x.^10;x.^11;x.^12;x.^13;x.^14;x.^15;x.^16;x.^17;x.^18;x.^19;x.^20;x.^21;x.^22;x.^23;x.^24;x.^25;x.^26;x.^27;x.^28;x.^29;x.^30;x.^31;x.^32;x.^33;x.^34;x.^35;x.^36;x.^37;x.^38;x.^39;x.^40;x.^41;x.^42;x.^43;x.^44;x.^45;x.^46;x.^47;x.^48;x.^49;x.^50;x.^51;x.^52;x.^53;x.^54;x.^55;x.^56;x.^57;x.^58;x.^59;x.^60;x.^61;x.^62;x.^63;x.^64;x.^65;x.^66;x.^67;x.^68;x.^69;x.^70;x.^71;x.^72;x.^73;x.^74;x.^75;x.^76;x.^77;x.^78;x.^79;x.^80;x.^81;x.^82;x.^83;x.^84;x.^85;x.^86;x.^87;x.^88;x.^89;x.^90;x.^91;x.^92;x.^93;x.^94;x.^95;x.^96;x.^97;x.^98;x.^99;x.^100]), 0, 1) - lam'*m1;

    m = 10;
    res_t = zeros(m, 6);
    res_a = zeros(m, 6);
    for i = 1:m
        
        lam = lam_his(389-i, :)';
        
        tic
        res1 = g(lam);
        t1 = toc;

        f_int = @(x, lam) exp(lam'*power(x, 0:100)');

        tic
        res2 = guassian_quadrature_100(@(x) f_int(x, lam) - lam'*m1);
        t2 = toc;

        tic
        res3 = guassian_quadrature_500(@(x) f_int(x, lam) - lam'*m1);
        t3 = toc;

        tic
        res4 = guassian_quadrature_1000(@(x) f_int(x, lam) - lam'*m1);
        t4 = toc;

        tic
        nn = 14;
        res5 = romberg_integral(@(x) f_int(x, lam) - lam'*m1, 0, 1, nn);
        t5 = toc;

        tic
        nn = 18;
        res6 = romberg_integral(@(x) f_int(x, lam) - lam'*m1, 0, 1, nn);
        t6 = toc;

        tic
        nn = 20;
        res7 = romberg_integral(@(x) f_int(x, lam) - lam'*m1, 0, 1, nn);
        t7 = toc;
        
        res_t(i, :) = [t2, t3, t4, t5, t6, t7];
        res_a(i, :) = [res2, res3, res4, res5, res6, res7] - res1;
    end
end

    function [res] = guassian_quadrature_100(f)
    % n is the number of nodes
    
        global w_100 x_100

        ff = arrayfun(f, x_100);
        res = 0.5 .* ff' * w_100;
    end


    function [res] = guassian_quadrature_500(f)
    % n is the number of nodes
    
        global w_500 x_500
        
        ff = arrayfun(f, x_500);
        res = 0.5 .* ff' * w_500;
    end

    function [res] = guassian_quadrature_1000(f)
    % n is the number of nodes
        
        global w_1000 x_1000
        
        ff = arrayfun(f, x_1000);
        res = 0.5.* ff' * w_1000;
    end

    
    function [res] = romberg_integral (f, start_, end_, n)
    % n is the size of Romberg table

        % initialization
        h = end_ - start_;
        R = zeros(2, n);
        R(1, 1) = h/2 * (f(start_) + f(end_));

        for i = 2:n
            
            sum_ = sum(arrayfun(f, start_ + ((1:2^(i-2))-0.5)*h));
            R(2, 1) = 1/2 * (R(1, 1) + h*sum_);

            for j = 2:i
                R(2, j) = R(2, j-1) + (R(2, j-1)-R(1, j-1)) / (4^(j-1)-1);
            end

            h = h/2;
            R(1, :) = R(2, :);
        end

        res = R(2, n);
    end
