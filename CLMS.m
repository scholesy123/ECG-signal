function [w_final, y_final, e_final, all_w] = CLMS(N0, miu1, miu2, miu_a, a_plus, w1_init, w2_init, a_init, x, d)
    order = length(w1_init);
    y1 = zeros(length(x), 1);
    y2 = zeros(length(x), 1);
    y = zeros(length(x), 1);
    y_final = zeros(length(x), 1);
    e1 = zeros(length(x), 1);
    e2 = zeros(length(x), 1);
    e = zeros(length(x), 1);
    e_final = zeros(length(x), 1);
    w1 = w1_init;
    w2 = w2_init;
    lambda = 1 / (1 + exp(-a_init));
    a = a_init;
    all_w = zeros(order, length(x) - order);
    for n = 1 : (length(x) - order)
        y1(n + order) = w1' * x(n : n + order - 1);
        y2(n + order) = w2' * x(n : n + order - 1);
        e1(n + order) = d(n + order) - y1(n + order);
        e2(n + order) = d(n + order) - y2(n + order);
        y(n + order) = lambda * y1(n + order) + (1 - lambda) * y2(n + order);
        e(n + order) = d(n + order) - y(n + order);
        w1 = w1 + miu1 * e1(n + order) * x(n : n + order - 1);
        a = a + miu_a * sign(e(n + order) * (e2(n + order) - e1(n + order)));
        lambda = 1 / (1 + exp(-a));
        if a < -a_plus
            a = -a_plus;
            lambda = 0;
        end
        if a > a_plus
            a = a_plus;
            lambda = 1;
            if (mod(n - 1, N0)) == 0
                w2 = w1;
            end
        else
            w2 = w2 + miu2 * e2(n + order) * x(n : n + order - 1);
        end
        w = lambda * w1 + (1 - lambda) * w2;
        all_w(:, n) = w;
    end
    w_final = w;
    for i = 1 : (length(x) - order)
        y_final(i + order) = w_final' * x(i : i + order - 1);
        e_final(i + order) = d(i + order) - y_final(i + order);
    end
end
        



