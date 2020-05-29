function sde(func_no, runs, bounds, optima, tolerance, popSize, Max_gen, Rs, t_dis)
    % function cde(runs,func_no)
    warning off;

    global xl xu
    F = 0.9;
    CR = 0.1;
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

    % [sortval,sortindex]=sort(val,'descend');
    % popsort=pop(sortindex,:);
    % valsort=val(sortindex);
    % spop.species=[];
    % i=1;
    %
    %
    %   while i<=size(popsort,1)-1,
    %     dist=zeros(size(popsort,1),1);
    %     dist(i+1:size(popsort,1),:)=sum((ones(size(popsort,1)-i,1)*popsort(i,:)-popsort(i+1:size(popsort,1),:)).^2,2)<Rs.^2;
    %     spop(i).pop=[popsort(i,:);popsort(dist==1,:)];
    %     spop(i).val=[val(i);valsort(dist==1)'];
    %     popsort=popsort(dist==0,:);
    %     valsort=valsort(dist==0);
    %     spop(i).species=popsort(i,:);
    %     spop(i).speciesval=valsort(i);
    %     i=i+1;
    %   end

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

        % sizepopsort=size(popsort,1)

        while i <= size(popsort, 1) - 1,
            dist = zeros(size(popsort, 1), 1);
            dist(i + 1:size(popsort, 1), :) = sum((ones(size(popsort, 1) - i, 1) * popsort(i, :) - popsort(i + 1:size(popsort, 1), :)).^2, 2) < Rs.^2;
            spop(i).pop = [popsort(i, :); popsort(dist == 1, :)];
            spop(i).val = [valsort(i); valsort(dist == 1)'];

            spop(i).species = popsort(i, :);
            spop(i).speciesval = valsort(i);
            popsort = popsort(dist == 0, :);
            valsort = valsort(dist == 0);
            i = i + 1;
        end

        if size(popsort, 1) - i + 1 ~= 0
            spop(i).species = popsort(i, :);
            spop(i).speciesval = valsort(i);
            spop(i).pop = spop(i).species;
            spop(i).val = spop(i).speciesval;
        end

        %   pop1=[];
        %   val1=[];
        %      for i=1:size(spop,2)
        %          pop1=[pop1;spop(i).pop];
        %          val1=[val1,spop(i).val'];
        %      end
        %
        %   sizepop1=size(pop1,1)

        for i = 1:size(spop, 2)

            while size(spop(i).pop, 1) < 10
                temprang = sqrt(Rs^2 / D);
                newpop = spop(i).species + temprang .* (-1 + 2 .* rand(1, D));
                newpop = (newpop < xl) .* (xl + rand(1, D) .* (xu - xl)) + (newpop >= xl) .* newpop;
                newpop = (newpop > xu) .* (xl + rand(1, D) .* (xu - xl)) + (newpop <= xu) .* newpop;
                newpopval = eobj(newpop, func_no);
                spop(i).pop = [spop(i).pop; newpop];
                spop(i).val = [spop(i).val; newpopval];
                nfeval = nfeval + 1;
            end

            for j = 1:size(spop(i).pop, 1)
                popold = spop(i).pop(j, :);
                newpop1 = spop(i).pop;
                newval1 = spop(i).val;
                bm = spop(i).species;
                ui(j, 1:D) = DE(popold, newpop1, bm, st, F, CR, D, size(newpop1, 1));
                tempval(j, 1) = eobj(ui(j, :), func_no);
                nfeval = nfeval + 1;

                if tempval(j, 1) == spop(i).speciesval
                    ui(j, :) = XRmin + (XRmax - XRmin) .* rand(1, D);
                    tempval(j, 1) = eobj(ui(j, :), func_no);
                    nfeval = nfeval + 1;

                end

                if tempval(j, 1) > spop(i).val(j)
                    spop(i).val(j) = tempval(j, 1);
                    spop(i).pop(j, :) = ui(j, :);
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

        %       traceInfo2(iter,1)=nfeval;
        %       traceInfo2(iter,2)=max(val);      % recording the best fitness
        %       traceInfo2(iter,3)=mean(val);     % recording the Avg fitness
        %       traceInfo2(iter,4)=std(val);      % recording the StdDev of fitness
        %       foundIn=zeros(numOpt,1);
        %       peaks=zeros(numOpt,1);
        %        distance=10000*ones(numOpt,1);
        %       for u=1:popSize,
        %
        %           checkdis=sum((ones(numOpt,1)*pop(u,:)-optima).^2,2);
        %
        %         distIn=checkdis<=tolerance.^2;
        % %         distance=sqrt(checkdis);
        %
        %         distance=(distance<sqrt(checkdis)).*distance+(distance>=sqrt(checkdis)).*sqrt(checkdis);
        %
        %         foundIn=max(foundIn,distIn);
        %         peaks=max(peaks,distIn.*val(u));
        %       end
        %       traceInfo2(iter,5)=sum(foundIn)/numOpt;               % recording the success rate
        % %       if sum(foundIn)==0
        %         traceInfo2(iter,6)=0;
        % %       else
        % %         traceInfo(iter,6)=sum(peaks)/sum(OptFit.*foundIn); % recording MPR
        % %       end
        %     traceInfo2(iter,7:numOpt+6)=distance;
        %       if traceInfo2(iter,5)==1
        %           break;
        %       end
        if traceInfo(iter, 5) == 1
            break;
        end

    end %---end while ((iter < Max_Gen) ...

    finalpeak = [peakslocation, peaks];
    endPop = [pop val'];
    dlmwrite(strcat('sde_info', char(num2str(func_no)), '_', char(num2str(runs)), '.txt'), traceInfo, 'newline', 'pc');
    dlmwrite(strcat('sde_result', char(num2str(func_no)), '_', char(num2str(runs)), '.txt'), endPop, 'newline', 'pc');
    dlmwrite(strcat('sde_peak', char(num2str(func_no)), '_', char(num2str(runs)), '.txt'), finalpeak, 'newline', 'pc');
    clear all
