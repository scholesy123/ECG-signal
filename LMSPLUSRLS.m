function [w_final, y_final, e_final, all_w] = LMSPLUSRLS(delta, miu, lambda, miu_a, a_plus, w1_init, w2_init, a_init, x, d)
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
    eta = 1 / (1 + exp(-a_init));
    a = a_init;
    P = 1 / delta * eye(order);
    all_w = zeros(order, length(x) - order);
    for n = 1 : (length(x) - order)
        y1(n + order) = w1' * x(n : n + order - 1);
        y2(n + order) = w2' * x(n : n + order - 1);
        e1(n + order) = d(n + order) - y1(n + order);
        e2(n + order) = d(n + order) - y2(n + order);
        y(n + order) = eta * y1(n + order) + (1 - eta) * y2(n + order);
        e(n + order) = d(n + order) - y(n + order);
        w1 = w1 + miu * e1(n + order) * x(n : n + order - 1);
        pai = 1 / lambda * P * x(n : n + order - 1);
        k = pai / (1 + x(n : n + order - 1)' * pai);
        w2 = w2 + k * e2(n + order);
        P = 1 / lambda * P - 1 / lambda * k * x(n : n + order - 1)' * P;
        e_alpha = e(n + order) * (y1(n + order) - y2(n + order));
        a = a + miu_a * e_alpha * eta * (1 - eta);
        eta = 1 / (1 + exp(-a));
        if a < -a_plus
            a = -a_plus;
            eta = 0;
        elseif a > a_plus
            a = a_plus;
            eta = 1;
        end
        w = eta * w1 + (1 - eta) * w2;
        all_w(:, n) = w;
    end
    w_final = w;
    for i = 1 : (length(x) - order)
        y_final(i + order) = w_final' * x(i : i + order - 1);
        e_final(i + order) = d(i + order) - y_final(i + order);
    end
end
        






