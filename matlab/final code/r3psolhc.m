% function [gbest,gbestval,fitcount,R1,R2]= PSO_func(fhd,Max_Gen,Max_FES,Particle_Number,Dimension,VRmin,VRmax,varargin)
function r2pso(func_no, runs, bounds, optima, tolerance, Particle_Number, Max_Gen, t_dis)

    Dimension = size(bounds, 1);
    Max_FES = Max_Gen * Particle_Number;
    VRmin = bounds(:, 1)';
    VRmax = bounds(:, 2)';
    numOpt = size(optima, 1);
    OptFit = zeros(numOpt, 1);
    OptFit = eobj(optima, func_no);

    % rand('state',sum(100*clock));
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
    nbest = pos;
    nbestval = e; %initialize the nbest and the nbest's fitness value
    i = 1;

    while i <= ps

        if i < ps - 1

            if pbestval(i + 1) > nbestval(i)
                tempval = pbestval(i + 1);
                temp = pbest(i + 1);
            else
                tempval = pbestval(i);
                temp = pbest(i);
            end

            if pbestval(i + 2) > tempval
                tempval = pbestval(i + 2);
                temp = pbest(i + 2);
            end

            nbestval(i) = tempval;
            nbestval(i + 1) = tempval;
            nbestval(i + 2) = tempval;
            nbest(i) = temp;
            nbest(i + 1) = temp;
            nbest(i + 2) = temp;
        elseif i == ps - 1

            if pbestval(i + 1) > nbestval(i)
                tempval = pbestval(i + 1);
                temp = pbest(i + 1);
            else
                tempval = pbestval(i);
                temp = pbest(i);
            end

            nbestval(i) = tempval;
            nbestval(i + 1) = tempval;
            nbest(i) = temp;
            nbest(i + 1) = temp;
        end

        i = i + 3;
    end

    [gbestval, gbestid] = max(pbestval);
    gbest = pbest(gbestid, :); %initialize the gbest and the gbest's fitness value
    gbestrep = repmat(gbest, ps, 1);
    g_res(1) = gbestval;

    % tmp1=abs(repmat(gbest,ps,1)-pos)+abs(pbest-pos);
    % temp1=ones(ps,1);
    % temp2=ones(ps,1);

    % for kkk=1:D
    %     temp1=temp1.*tmp1(:,kkk);
    %     temp2=temp2.*(max(pbest(:,kkk))-min(pbest(:,kkk)));
    % end
    % R1(1)=mean(temp1);
    % R2(1)=mean(temp2);

    for i = 2:me

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
            %     if pbestval(k)>gbestval
            %     gbest=pbest(k,:);
            %     gbestval=pbestval(k);
            %     gbestrep=repmat(gbest,ps,1);%update the gbest
            %     end

        end

        nbest = pbest;
        nbestval = pbestval;

        j = 1;

        while j <= ps

            if j < ps - 1

                if pbestval(j + 1) > nbestval(j)
                    tempval = pbestval(j + 1);
                    temp = pbest(j + 1);
                else
                    tempval = pbestval(j);
                    temp = pbest(j);
                end

                if pbestval(j + 2) > tempval
                    tempval = pbestval(j + 2);
                    temp = pbest(j + 2);
                end

                nbestval(j) = tempval;
                nbestval(j + 1) = tempval;
                nbestval(j + 2) = tempval;
                nbest(j) = temp;
                nbest(j + 1) = temp;
                nbest(j + 2) = temp;
            elseif j == ps - 1

                if pbestval(j + 1) > nbestval(j)
                    tempval = pbestval(j + 1);
                    temp = pbest(j + 1);
                else
                    tempval = pbestval(j);
                    temp = pbest(j);
                end

                nbestval(j) = tempval;
                nbestval(j + 1) = tempval;
                nbest(j) = temp;
                nbest(j + 1) = temp;
            end

            j = j + 3;
        end

        i

        traceInfo(i, 1) = i * ps;
        traceInfo(i, 2) = max(pbestval); % recording the best fitness
        traceInfo(i, 3) = mean(pbestval); % recording the Avg fitness
        traceInfo(i, 4) = std(pbestval); % recording the StdDev of fitness

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

        traceInfo(i, 5) = foundIn / numOpt;
        traceInfo(i, 6) = 0;
        traceInfo(i, 7:numOpt + 6) = peaks;

        %       traceInfo(i,1)=i*ps;
        %       traceInfo(i,2)=max(pbestval);      % recording the best fitness
        %       traceInfo(i,3)=mean(pbestval);     % recording the Avg fitness
        %       traceInfo(i,4)=std(pbestval);      % recording the StdDev of fitness
        %       foundIn=zeros(numOpt,1);
        %       peaks=zeros(numOpt,1);
        %        distance=10000*ones(numOpt,1);
        %       for u=1:ps
        %
        %           checkdis=sum((ones(numOpt,1)*pbest(u,:)-optima).^2,2);
        %
        %         distIn=checkdis<=tolerance.^2;
        % %         distance=sqrt(checkdis);
        %
        %         distance=(distance<sqrt(checkdis)).*distance+(distance>=sqrt(checkdis)).*sqrt(checkdis);
        %
        %         foundIn=max(foundIn,distIn);
        % %         peaks=max(peaks,distIn.*pbestval(u));
        %       end
        %       traceInfo(i,5)=sum(foundIn)/numOpt;               % recording the success rate
        % %       if sum(foundIn)==0
        %         traceInfo(i,6)=0;
        % %       else
        % %         traceInfo(iter,6)=sum(peaks)/sum(OptFit.*foundIn); % recording MPR
        % %       end
        %     traceInfo(i,7:numOpt+6)=distance;
        if traceInfo(i, 5) == 1;
            break;
        end

        if fitcount >= Max_FES
            break;
        end

    end

    finalpeak = [peakslocation, peaks];
    fresult = [pbest, pbestval];
    dlmwrite(strcat('r3psolhc_info', char(num2str(func_no)), '_', char(num2str(runs)), '.txt'), traceInfo, 'newline', 'pc');
    dlmwrite(strcat('r3psolhc_result', char(num2str(func_no)), '_', char(num2str(runs)), '.txt'), fresult, 'newline', 'pc');
    dlmwrite(strcat('r3psolhc_peak', char(num2str(func_no)), '_', char(num2str(runs)), '.txt'), finalpeak, 'newline', 'pc');
    clear all
