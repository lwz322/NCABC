function cde(func_no, runs, bounds, optima, tolerance, popSize, Max_gen, t_dis)
    % function cde(runs,func_no)
    warning off;

    global xl xu
    F = 0.9;
    CR = 0.2;
    st = 2;

    % NP=popSize/10;
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

        for j = 1:NP

            popold = pop(j, :);
            newpop1 = pop;
            newval1 = val;
            bm = DE_gbest;
            ui(j, 1:D) = DE(popold, newpop1, bm, st, F, CR, D, size(newpop1, 1));

            nfeval = nfeval + 1;
        end

        %-----Select which vectors are allowed to enter the new population------------

        tempval = eobj(ui, func_no);

        for i = 1:NP

            [temp j] = min(sqrt(sum((ones(NP, 1) * ui(i, :) - pop).^2, 2)));

            if val(j) < tempval(i)
                pop(j, :) = ui(i, :);
                val(j) = tempval(i);

                if (tempval(i) > DE_gbestval)% if competitor better than the best one ever
                    DE_gbestval = tempval(i); % new best value
                    DE_gbest = ui(i, :); % new best parameter vector ever

                end

            end

        end %---end for imember=1:NP

        DE_gbestval
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

            if abs(val(minindex) - OptFit(u)) <= tolerance & minval <= t_dis^2

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

    end %---end while ((iter < Max_Gen) ...

    finalpeak = [peakslocation, peaks];
    endPop = [pop val'];
    dlmwrite(strcat('cde_info', char(num2str(func_no)), '_', char(num2str(runs)), '.txt'), traceInfo, 'newline', 'pc');
    dlmwrite(strcat('cde_result', char(num2str(func_no)), '_', char(num2str(runs)), '.txt'), endPop, 'newline', 'pc');
    dlmwrite(strcat('cde_peak', char(num2str(func_no)), '_', char(num2str(runs)), '.txt'), finalpeak, 'newline', 'pc');
    clear all
