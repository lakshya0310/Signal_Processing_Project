

Ts = 0.001;
N  = 1000;
Ks = 1:4;
R  = 50;

signals = {
    @(t) cos(2*pi*250*t);
    @(t) 0.6*cos(2*pi*40*t) + 0.4*cos(2*pi*90*t);
    @(t) sinc(100*(t-0.3));
};
sig_names = {
    'Signal 1: cos(2\pi 50 t)';
    'Signal 2: 0.6cos(2\pi40t)+0.4cos(2\pi90t)';
    'Signal 3: sinc(100(t-0.3))';
};

for s = 1:numel(signals)
    x_fun = signals{s};
    fprintf('\n%s\n', sig_names{s});
    fprintf("Scenario (a): MSE vs K (spline reconstruction)\n");
    fprintf("============================================\n");
    fprintf("  K    |   mean MSE over %d runs\n", R);
    fprintf("--------+------------------------\n");

    MSEs = zeros(size(Ks));

    t_cont = linspace(0, (N-1)*Ts, 5000);
    x_cont = x_fun(t_cont);

    figure('Name', sig_names{s});

    for idx = 1:numel(Ks)
        K = Ks(idx);
        mse_acc = 0;
        x_rec_last = [];
        t_ideal_last = [];

        for r = 1:R
            [x_jit, k, t_jit, t_ideal] = jittered_sampling(x_fun, Ts, N, K);
            x_rec  = jittered_spline_reconstruct(t_ideal, t_jit, x_jit);
            x_true = x_fun(t_ideal).';

            mse_acc = mse_acc + mean((x_true - x_rec).^2);

            if r == R
                x_rec_last   = x_rec;
                t_ideal_last = t_ideal;
            end
        end

        MSEs(idx) = mse_acc / R;
        fprintf("  %d    |   %.20f\n", K, MSEs(idx));

        subplot(4,1,idx);
        plot(t_cont, x_cont, 'LineWidth', 1.1); hold on;
        plot(t_ideal_last, x_rec_last, '--', 'LineWidth', 1.1);
        grid on;
        if idx == 1
            title(sprintf('%s  (K = %d, MSE = %.20e)', sig_names{s}, K, MSEs(idx)));
        else
            title(sprintf('K = %d, MSE = %.20e', K, MSEs(idx)));
        end
        ylabel('Amplitude');
        if idx == numel(Ks)
            xlabel('t (s)');
        end
        legend('True continuous x(t)', 'Reconstructed', 'Location', 'best');
    end

    fprintf("============================================\n");
end


function [x_hat, k, t_jittered, t_ideal] = jittered_sampling(x_fun, Ts, N, K)
    Delta   = Ts / 10;
    n       = 0:N-1;
    t_ideal = n * Ts;
    k = randi([-K, K], 1, N);
    t_jittered = t_ideal + k * Delta;
    x_hat = x_fun(t_jittered);
end

function x_rec = jittered_spline_reconstruct(t_eval, t_jittered, x_hat)
    t_eval     = t_eval(:);
    t_jittered = t_jittered(:).';
    x_hat      = x_hat(:).';
    x_rec = spline(t_jittered, x_hat, t_eval);
end