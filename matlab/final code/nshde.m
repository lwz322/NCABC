function nshde(func_no, runs, bounds, optima, tolerance, popSize, Max_gen, Rs, t_dis)
    % function cde(runs,func_no)
    warning off;

    global xl xu
    F = 0.9;
    CR = 0.1;
    st = 2;

    % NP=popSize/10;
    NP = popSize;

    Msize = 2 * NP;

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

    for i = 1:Msize
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

    for i = 2:Msize% check the remaining members

        val(i) = eobj(pop(i, :), func_no);

        nfeval = nfeval + 1;

        if (val(i) > DE_gbestval)% if member is better
            ibest = i; % save its location
            DE_gbestval = val(i);
        end

    end

    DE_gbest = pop(ibest, :); % best member of current iteration
    bestvalit = DE_gbestval; % best value of current iteration

    mpop = pop;
    mval = val';

    iter = 0;

    while nfeval < Max_FES

        %     [sortval,sortindex]=sort(mval,'descend');
        %     popsort=mpop(sortindex,:);
        %     valsort=mval(sortindex);

        spop = sharing(mpop, mval, Rs);
        mpop = [];
        mval = [];
        pval = [];

        for i = 1:size(spop, 2)
            mpop = [mpop; spop(i).pop];
            mval = [mval; spop(i).val];
            pval = [pval; spop(i).newval];
        end

        [sortpval, sortpindex] = sort(pval, 'descend');
        valsort = mval(sortpindex);
        popsort = mpop(sortpindex, :);
        pvalsort = pval(sortpindex);

        if size(popsort, 1) > Msize
            popsort = popsort(1:Msize, :);
            valsort = valsort(1:Msize);
            pvalsort = pvalsort(1:Msize);
        end

        for j = 1:NP

            if NP <= 200
                numb = 5 + 5 * ((Max_FES - nfeval) / Max_FES);
            else
                numb = 20 + 30 * ((Max_FES - nfeval) / Max_FES);
            end

            %         numb=NP/10;

            bm = popsort(j, :);

            [temp k] = sort(sqrt(sum((ones(size(popsort, 1), 1) * bm - popsort).^2, 2)));

            popold = popsort(j, :);

            newpop1 = popsort(k(1:numb), :);

            ui(j, 1:D) = DE(popold, newpop1, bm, st, F, CR, D, size(newpop1, 1));

            nfeval = nfeval + 1;
        end

        tempval = eobj(ui, func_no);

        for i = 1:NP

            pp = randperm(Msize);
            qq = popsort(pp(1:Msize), :);

            [temp j] = min(sqrt(sum((ones(Msize, 1) * ui(i, :) - qq).^2, 2)));

            if valsort(pp(j)) < tempval(i)
                popsort(pp(j), :) = ui(i, :);
                valsort(pp(j)) = tempval(i);

            end

        end

        mpop = popsort;
        mval = valsort;

        max(mval)

        iter = iter + 1

        traceInfo(iter, 1) = nfeval;
        traceInfo(iter, 2) = max(mval); % recording the best fitness
        traceInfo(iter, 3) = mean(mval); % recording the Avg fitness
        traceInfo(iter, 4) = std(mval); % recording the StdDev of fitness

        peaks = -10000 * ones(numOpt, 1);
        distance = 10000 * ones(numOpt, 1);
        peakslocation = 10000 * ones(numOpt, D);

        for u = 1:numOpt
            checkdis = sum((ones(size(mpop, 1), 1) * optima(u, :) - mpop).^2, 2);
            [minval, minindex] = min(checkdis);

            if abs(mval(minindex) - OptFit(u)) <= tolerance & minval <= t_dis

                if mval(minindex) > peaks(u)
                    peaks(u) = mval(minindex);
                    peakslocation(u, :) = mpop(minindex, :);
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
    endPop = [mpop mval];
    dlmwrite(strcat('nshde_info', char(num2str(func_no)), '_', char(num2str(runs)), '.txt'), traceInfo, 'newline', 'pc');
    dlmwrite(strcat('nshde_result', char(num2str(func_no)), '_', char(num2str(runs)), '.txt'), endPop, 'newline', 'pc');
    dlmwrite(strcat('nshde_peak', char(num2str(func_no)), '_', char(num2str(runs)), '.txt'), finalpeak, 'newline', 'pc');
    clear all
