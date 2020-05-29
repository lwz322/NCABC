% function [gbest,gbestval,fitcount,R1,R2]= PSO_func(fhd,Max_Gen,Max_FES,Particle_Number,Dimension,VRmin,VRmax,varargin)
function ferpso(func_no, runs, bounds, optima, tolerance, Particle_Number, Max_Gen, t_dis)

    Dimension = size(bounds, 1);
    Max_FES = Max_Gen * Particle_Number;
    VRmin = bounds(:, 1)';
    VRmax = bounds(:, 2)';
    numOpt = size(optima, 1);
    OptFit = zeros(numOpt, 1);

    OptFit = eobj(optima, func_no);

    me = Max_Gen;
    ps = Particle_Number;
    D = Dimension;
    cc = [2.05 2.05]; %acceleration constants
    iwt = 0.729843788;

    if length(VRmin) == 1
        VRmin = repmat(VRmin, 1, D);
        VRmax = repmat(VRmax, 1, D);
    end

    mv = 0.5 * (VRmax - VRmin);
    VRmin = repmat(VRmin, ps, 1);
    VRmax = repmat(VRmax, ps, 1);
    Vmin = repmat(-mv, ps, 1);
    Vmax = -Vmin;
    pos = VRmin + (VRmax - VRmin) .* rand(ps, D);

    for i = 1:ps;
        e(i, 1) = eobj(pos(i, :), func_no);
    end

    fitcount = ps;
    vel = Vmin + 2 .* Vmax .* rand(ps, D); %initialize the velocity of the particles
    pbest = pos;
    pbestval = e; %initialize the pbest and the pbest's fitness value

    sos = sqrt(sum((bounds(:, 2)' - bounds(:, 1)').^2)); % get the size ofthe search space

    for gen = 2:me

        [gbestval, gbestid] = max(pbestval); % get the current gobal best value
        [gworstval, gworstid] = min(pbestval); % get the current gobal worst value
        sf = sos / (gbestval - gworstval); % calculate the current scale factor
        nbest = pbest; %inital the nbest position
        nbestval = pbestval; %initalthe nbest value

        % update the nbest and nbestval
        for i = 1:ps

            for j = 1:ps
                eucdis = sqrt(sum((pbest(i, :) - pbest(j, :)).^2, 2));

                if j == 1
                    fertmp = -inf;
                end

                if eucdis ~= 0
                    FER = sf * ((pbestval(j) - pbestval(i)) / eucdis);

                    if FER > fertmp
                        fertmp = FER;
                        nbest(i, :) = pbest(j, :);
                        nbestval(i) = pbestval(j);
                    end

                end

            end

        end

        for k = 1:ps

            aa(k, :) = cc(1) .* rand(1, D) .* (pbest(k, :) - pos(k, :)) + cc(2) .* rand(1, D) .* (nbest(k, :) - pos(k, :));
            vel(k, :) = iwt .* vel(k, :) + aa(k, :);
            vel(k, :) = (vel(k, :) > mv) .* mv + (vel(k, :) <= mv) .* vel(k, :);
            vel(k, :) = (vel(k, :) < (-mv)) .* (-mv) + (vel(k, :) >= (-mv)) .* vel(k, :);
            pos(k, :) = pos(k, :) + vel(k, :);
            pos(k, :) = ((pos(k, :) >= VRmin(1, :)) & (pos(k, :) <= VRmax(1, :))) .* pos(k, :) ...
                +(pos(k, :) < VRmin(1, :)) .* (VRmin(1, :) + 0.25 .* (VRmax(1, :) - VRmin(1, :)) .* rand(1, D)) + (pos(k, :) > VRmax(1, :)) .* (VRmax(1, :) - 0.25 .* (VRmax(1, :) - VRmin(1, :)) .* rand(1, D));

            e(k, 1) = eobj(pos(k, :), func_no);
            fitcount = fitcount + 1;
            tmp = (pbestval(k) > e(k));
            temp = repmat(tmp, 1, D);
            pbest(k, :) = temp .* pbest(k, :) + (1 - temp) .* pos(k, :);
            pbestval(k) = tmp .* pbestval(k) + (1 - tmp) .* e(k); %update the pbest

        end

        gen

        traceInfo(gen, 1) = gen * ps;
        traceInfo(gen, 2) = max(pbestval); % recording the best fitness
        traceInfo(gen, 3) = mean(pbestval); % recording the Avg fitness
        traceInfo(gen, 4) = std(pbestval); % recording the StdDev of fitness

        peaks = -10000 * ones(numOpt, 1);
        distance = 10000 * ones(numOpt, 1);
        peakslocation = 10000 * ones(numOpt, D);

        for u = 1:numOpt
            checkdis = sum((ones(ps, 1) * optima(u, :) - pbest).^2, 2);
            [minval, minindex] = min(checkdis);

            if abs(pbestval(minindex) - OptFit(u)) <= tolerance & minval <= t_dis

                if pbestval(minindex) > peaks(u)
                    peaks(u) = pbestval(minindex);
                    peakslocation(u, :) = pbest(minindex, :);
                end

            end

        end

        foundIn = 0;

        for u = 1:numOpt

            if peaks(u) >- 10000
                foundIn = foundIn + 1;
            end

        end

        traceInfo(gen, 5) = foundIn / numOpt;
        traceInfo(gen, 6) = 0;
        traceInfo(gen, 7:numOpt + 6) = peaks;

        if traceInfo(gen, 5) == 1;
            break;
        end

        if fitcount >= Max_FES
            break;
        end

    end

    finalpeak = [peakslocation, peaks];
    fresult = [pbest, pbestval];
    dlmwrite(strcat('ferpso_info', char(num2str(func_no)), '_', char(num2str(runs)), '.txt'), traceInfo, 'newline', 'pc');
    dlmwrite(strcat('ferpso_result', char(num2str(func_no)), '_', char(num2str(runs)), '.txt'), fresult, 'newline', 'pc');
    dlmwrite(strcat('ferpso_peak', char(num2str(func_no)), '_', char(num2str(runs)), '.txt'), finalpeak, 'newline', 'pc');
    clear all
