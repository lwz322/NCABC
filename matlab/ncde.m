% function [DE_gbest,DE_gbestval,DE_fitcount,DE_fit_cut,DE_get_flag]  = dealgorithmorigincompe(fname,VTR,Max_FES,D,XRmin,XRmax,Lbound,Ubound,NP,Max_Gen,F,CR,strategy,fun);
% function [endPop,traceInfo] = cde(func_no)
% function rts(runs,func_no)

function ncde(func_no, runs, bounds, optima, tolerance, popSize, Max_gen, t_dis)
    warning off;

    global xl xu
    F = 0.9;
    CR = 0.1;
    st = 2;

    warning off;

    NP = popSize;

    XRmin = bounds(:, 1)';
    XRmax = bounds(:, 2)';
    D = size(bounds, 1);
    Lbound = bounds(:, 1)';
    Ubound = bounds(:, 2)';

    Max_FES = NP * Max_gen;

    xl = bounds(:, 1)';
    xu = bounds(:, 2)';

    numOpt = size(optima, 1);
    OptFit = zeros(numOpt, 1);
    OptFit = eobj(optima, func_no);

    for i = 1:popSize
        pop(i, :) = XRmin + (XRmax - XRmin) .* rand(1, D);
    end

    val = zeros(1, popSize); % create and reset the "cost array"
    DE_gbest = zeros(1, D); % best population member ever
    nfeval = 0; % number of function evaluations

    %------Evaluate the best member after initialization----------------------

    ibest = 1; % start with first population member

    val(1) = eobj(pop(ibest, :), func_no);

    DE_gbestval = val(1); % best objective function value so far
    nfeval = nfeval + 1;

    for i = 2:popSize% check the remaining members

        val(i) = eobj(pop(i, :), func_no);

        nfeval = nfeval + 1;

        if (val(i) > DE_gbestval)% if member is better
            ibest = i; % save its location
            DE_gbestval = val(i);
        end

    end

    DE_gbest = pop(ibest, :); % best member of current iteration
    bestvalit = DE_gbestval; % best value of current iteration

    %------DE-Minimization---------------------------------------------
    %------popold is the population which has to compete. It is--------
    %------static through one iteration. pop is the newly--------------
    %------emerging population.----------------------------------------

    iter = 0;

    while nfeval < Max_FES

        %     popold = pop(k(1:NP),:);                   % save the old population
        for j = 1:NP
            %         numb=50;

            if NP <= 200
                numb = 5 + 5 * ((Max_FES - nfeval) / Max_FES);
            else
                numb = 20 + 30 * ((Max_FES - nfeval) / Max_FES);
            end

            %         numb=NP;

            bm = pop(j, :);
            %      ranindex=randperm(popSize);
            %
            %     bm=pop(ranindex(1),:);
            %
            [temp k] = sort(sqrt(sum((ones(popSize, 1) * bm - pop).^2, 2)));

            popold = pop(j, :);

            %         popold=bm;

            newpop1 = pop(k(1:numb), :);

            ui(j, 1:D) = DE(popold, newpop1, bm, st, F, CR, D, size(newpop1, 1));

            nfeval = nfeval + 1;
        end

        tempval = eobj(ui, func_no);

        for i = 1:NP

            pp = randperm(popSize);
            qq = pop(pp(1:NP), :);

            [temp j] = min(sqrt(sum((ones(NP, 1) * ui(i, :) - qq).^2, 2)));

            if val(pp(j)) < tempval(i)
                pop(pp(j), :) = ui(i, :);
                val(pp(j)) = tempval(i);

                if (tempval(i) > DE_gbestval)% if competitor better than the best one ever
                    DE_gbestval = tempval(i); % new best value
                    DE_gbest = ui(i, :); % new best parameter vector ever

                end

            end

        end %---end for imember=1:NP

        DE_gbestval;
        iter = iter + 1

        traceInfo(iter, 1) = iter * NP;
        traceInfo(iter, 2) = max(val); % recording the best fitness
        traceInfo(iter, 3) = mean(val); % recording the Avg fitness
        traceInfo(iter, 4) = std(val); % recording the StdDev of fitness

        peaks = -10000 * ones(numOpt, 1);
        distance = 10000 * ones(numOpt, 1);
        peakslocation = 10000 * ones(numOpt, D);

        for u = 1:numOpt
            checkdis = sum((ones(NP, 1) * optima(u, :) - pop).^2, 2);
            [minval, minindex] = min(checkdis);

            if abs(val(minindex) - OptFit(u)) <= tolerance & minval <= t_dis

                if val(minindex) > peaks(u)
                    peaks(u) = val(minindex);
                    peakslocation(u, :) = pop(minindex, :);
                end

            end

        end

        foundIn = 0;

        for u = 1:numOpt

            if peaks(u) >- 10000
                foundIn = foundIn + 1;
            end

        end

        traceInfo(iter, 5) = foundIn / numOpt;
        traceInfo(iter, 6) = 0;
        traceInfo(iter, 7:numOpt + 6) = peaks;

        if traceInfo(iter, 5) == 1;
            break;
        end

    end

    finalpeak = [peakslocation, peaks];
    endPop = [pop val'];
    dlmwrite(strcat('ncde_info', char(num2str(func_no)), '_', char(num2str(runs)), '.txt'), traceInfo, 'newline', 'pc');
    dlmwrite(strcat('ncde_result', char(num2str(func_no)), '_', char(num2str(runs)), '.txt'), endPop, 'newline', 'pc');
    dlmwrite(strcat('ncde_peak', char(num2str(func_no)), '_', char(num2str(runs)), '.txt'), finalpeak, 'newline', 'pc');
    clear all
