function ui = DE(popold, pop, bm, st, F, CR, n, NP)
    global xl xu
    r1 = round(rand * NP); r2 = round(rand * NP); r3 = round(rand * NP); r4 = round(rand * NP); r5 = round(rand * NP);
    while (r1 == j | r1 == 0), r1 = ceil(rand * NP); end
    while (r2 == j | r2 == r1 | r2 == 0), r2 = ceil(rand * NP); end
    while (r3 == j | r3 == r1 | r3 == r2 | r3 == 0), r3 = ceil(rand * NP); end
    while (r4 == j | r4 == r1 | r4 == r2 | r4 == r3 | r4 == 0), r4 = ceil(rand * NP); end
    while (r5 == j | r5 == r1 | r5 == r2 | r5 == r3 | r5 == r4 | r5 == 0), r5 = ceil(rand * NP); end
    pm1 = pop(r1, 1:n);
    pm2 = pop(r2, 1:n);
    pm3 = pop(r3, 1:n);
    pm4 = pop(r4, 1:n);
    pm5 = pop(r5, 1:n);
    rotd = (0:1:n - 1);

    mui = rand(1, n) < CR; % all random numbers < CR are 1, 0 otherwise
    if mui == zeros(1, n), nn = randperm(n); mui(nn(1)) = 1; end

    if st > 5
        st = st - 5;
        mui = sort(mui');
        nn = floor(rand .* n);

        if nn > 0
            rtd = rem(rotd + nn, n);
            mui(:) = mui(rtd + 1); %rotate column i by n
        end

        mui = mui';
    end

    mpo = mui < 0.5; % inverse mask to mui

    if (st == 1)% DE/best/1   6
        ui = bm + F * (pm1 - pm2); % differential variation
        ui = popold .* mpo + ui .* mui; % binomial crossover
    elseif (st == 2)% DE/rand/1   7
        ui = pm3 + F * (pm1 - pm2); % differential variation
        ui = popold .* mpo + ui .* mui; % crossover
    elseif (st == 3)% DE/rand-to-best/1    8
        ui = popold + F * (bm - popold) + F * (pm1 - pm2);
        ui = popold .* mpo + ui .* mui; % crossover
    elseif (st == 4)% DE/best/2           9
        ui = bm + F * (pm1 - pm2 + pm3 - pm4); % differential variation
        ui = popold .* mpo + ui .* mui; % crossover
    elseif (st == 5)% DE/rand/2           10
        ui = pm5 + F * (pm1 - pm2 + pm3 - pm4); % differential variation
        ui = popold .* mpo + ui .* mui; % crossover
    end

    if rand > 0.5
        ui = (ui < xl) .* xl + (ui >= xl) .* ui;
        ui = (ui > xu) .* xu + (ui <= xu) .* ui;
    else
        ui = (ui < xl) .* (xl + rand(1, n) .* (xu - xl)) + (ui >= xl) .* ui;
        ui = (ui > xu) .* (xl + rand(1, n) .* (xu - xl)) + (ui <= xu) .* ui;
    end
