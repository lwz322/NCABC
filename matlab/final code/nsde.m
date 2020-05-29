function sde(func_no, runs, bounds, optima, tolerance, popSize, Max_gen, Rs, t_dis)
    % function cde(runs,func_no)
    warning off;

    global xl xu
    F = 0.9;
    CR = 0.1;
    st = 2;

    % Rs=[];

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

    %------popold is the population which has to compete. It is--------
    %------static through one iteration. pop is the newly--------------
    %------emerging population.----------------------------------------

    iter = 0;

    while nfeval < Max_FES

        [sortval, sortindex] = sort(val, 'descend');
        popsort = pop(sortindex, :);
        valsort = val(sortindex);
        clear spop;
        i = 1;

        if size(popsort, 1) > NP
            popsort = popsort(1:NP, :);
            valsort = valsort(1:NP);
        end

        for i = 1:(NP / 5)

            [temp k] = sort(sqrt(sum((ones(size(popsort, 1), 1) * popsort(1, :) - popsort).^2, 2)));
            spop(i).species = popsort(1, :);
            spop(i).speciesval = valsort(1);
            checker = ones(size(popsort, 1), 1);
            checker(k(1:5), :) = 0;
            spop(i).pop = popsort(checker == 0, :);
            spop(i).val = valsort(checker == 0)';
            popsort = popsort(checker == 1, :);
            valsort = valsort(checker == 1);
        end

        for i = 1:size(spop, 2)

            for j = 1:size(spop(i).pop, 1)
                popold = spop(i).pop(j, :);
                newpop1 = spop(i).pop;
                newval1 = spop(i).val;
                bm = spop(i).species;
                ui(j, 1:D) = DE(popold, newpop1, bm, st, F, CR, D, size(newpop1, 1));
                tempval(j, 1) = eobj(ui(j, :), func_no);
                nfeval = nfeval + 1;

                checkdis = sum((ones(size(spop(i).pop, 1), 1) * ui(j, :) - spop(i).pop).^2, 2);
                [minval, minindex] = min(checkdis);

                if tempval(j, 1) > spop(i).val(minindex)
                    spop(i).val(minindex) = tempval(j, 1);
                    spop(i).pop(minindex, :) = ui(j, :);
                end

            end

        end

        pop = [];
        val = [];

        for i = 1:size(spop, 2)
            pop = [pop; spop(i).pop];
            val = [val, spop(i).val'];
        end

        iter = iter + 1

        traceInfo(iter, 1) = nfeval;
        traceInfo(iter, 2) = max(val); % recording the best fitness
        traceInfo(iter, 3) = mean(val); % recording the Avg fitness
        traceInfo(iter, 4) = std(val); % recording the StdDev of fitness

        peaks = -10000 * ones(numOpt, 1);
        distance = 10000 * ones(numOpt, 1);
        peakslocation = 10000 * ones(numOpt, D);

        for u = 1:numOpt
            checkdis = sum((ones(size(pop, 1), 1) * optima(u, :) - pop).^2, 2);
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

        if traceInfo(iter, 5) == 1
            break;
        end

    end %---end while ((iter < Max_Gen) ...

    finalpeak = [peakslocation, peaks];
    endPop = [pop val'];
    dlmwrite(strcat('nsde_info', char(num2str(func_no)), '_', char(num2str(runs)), '.txt'), traceInfo, 'newline', 'pc');
    dlmwrite(strcat('nsde_result', char(num2str(func_no)), '_', char(num2str(runs)), '.txt'), endPop, 'newline', 'pc');
    dlmwrite(strcat('nsde_peak', char(num2str(func_no)), '_', char(num2str(runs)), '.txt'), finalpeak, 'newline', 'pc');
    clear all
