

Tstep = 0.001;
L      = 1000;
p_show = 0.05;

sigSet = {
    @(t) cos(2*pi*50*t);
    @(t) 0.6*cos(2*pi*40*t) + 0.4*cos(2*pi*90*t);
    @(t) sinc(100*(t-0.3));
};

sigNames = {
    'Signal 1: cos(2\pi50t)';
    'Signal 2: 0.6cos(2\pi40t)+0.4cos(2\pi90t)';
    'Signal 3: sinc(100(t-0.3))';
};

for s = 1:length(sigSet)
    g = sigSet{s};

    [obsSeq, idxVec, trueSeq] = gen_missing(g, Tstep, L, p_show);

    [estSeq, tGrid, yRec] = interp_missing(obsSeq, Tstep);

    yTrueCont = g(tGrid);

    mse_demo = mean((trueSeq - estSeq).^2);
    fprintf('%s\n', sigNames{s});
    fprintf('Demo p = %.2f, MSE = %.6e\n\n', p_show, mse_demo);

    figure;
    plot(tGrid, yTrueCont, 'b', 'LineWidth', 1.2); hold on;
    plot(tGrid, yRec, 'r--', 'LineWidth', 1.2);
    title(sprintf('%s  |  p = %.2f, MSE = %.3e', ...
        sigNames{s}, p_show, mse_demo));
    xlabel('t'); ylabel('Amplitude');
    legend('True','Estimated');
    grid on;
end


pVals = 0.01:0.01:0.10;
Rtrials = 50;
errMat = zeros(length(sigSet), length(pVals));

for s = 1:length(sigSet)
    g = sigSet{s};

    fprintf('\n%s\n', sigNames{s});
    fprintf('Scenario (b) MSE vs p\n');
    fprintf('=====================================\n');
    fprintf('   p    |   mean MSE (%d trials)\n', Rtrials);
    fprintf('--------+------------------------------\n');

    for pi = 1:length(pVals)
        pp = pVals(pi);
        acc = 0;

        for r = 1:Rtrials
            [obsSeq, idxVec, trueSeq] = gen_missing(g, Tstep, L, pp);
            estOnly = interp_missing(obsSeq, Tstep);
            acc = acc + mean((trueSeq - estOnly).^2);
        end

        errMat(s, pi) = acc / Rtrials;
        fprintf(' %.2f   |   %.6e\n', pp, errMat(s,pi));
    end

    fprintf('=====================================\n');
end

figure; hold on;
fmt = {'-o','-s','-^'};
for s = 1:length(sigSet)
    plot(pVals, errMat(s,:), fmt{s}, 'LineWidth', 1.2);
end
xlabel('Missing probability p');
ylabel('MSE');
yscale log
title('Scenario (b): MSE vs p');
legend(sigNames, 'Location', 'northwest');
grid on;
hold off;


function [xmiss, idx, xideal] = gen_missing(funh, Ts, N, p)
    idx   = (0:N-1).';
    tvals = idx*Ts;
    xideal = funh(tvals);
    mask  = rand(N,1) > p;
    xmiss = xideal;
    xmiss(~mask) = NaN;
end

function [xest, tC, yC] = interp_missing(xin, Ts)
    xin = xin(:);
    L   = length(xin);
    xfilled = fillmissing(xin,'linear','EndValues','nearest');

    if nargout == 1
        xest = xfilled;
        return;
    end

    n = 0:L-1;
    tS = n*Ts;
    tC = linspace(tS(1), tS(end), 5000);

    xrow = xfilled(:).';
    trow = tS(:).';

    yC = zeros(length(tC),1);
    for m = 1:length(tC)
        tau = (tC(m) - trow) / Ts;
        yC(m) = sum(xrow .* sinc(tau));
    end

    xest = xfilled;
end
