function simulation(test_func_num, niching_method_num, runs)
    % simulation control function

    errorflag = 0;

    n = nargin;

    if n < 2
        disp('Insufficient arguements')
    end

    if test_func_num == 1
        bounds = [0 20];
        optima = 20;
        tolerance = 0.05;
        NP = 50;
        Max_Gen = 200;
        Rs = 0.5;
        t_dis = 0.1;
    elseif test_func_num == 2
        bounds = [0 20];
        optima = 20;
        tolerance = 0.05;
        NP = 50;
        Max_Gen = 200;
        Rs = 0.5;
        t_dis = 0.1;
    elseif test_func_num == 3
        bounds = [0 30];
        optima = [0 30]';
        tolerance = 0.05;
        NP = 50;
        Max_Gen = 200;
        Rs = 0.5;
        t_dis = 0.1;
    elseif test_func_num == 4
        bounds = [0 1];
        optima = [0.1 0.3 0.5 0.7 0.9]';
        tolerance = 0.000001;
        NP = 50;
        Max_Gen = 200;
        Rs = 0.01;
        t_dis = 0.01;
    elseif test_func_num == 5
        bounds = [0 1];
        optima = 0.1;
        tolerance = 0.000001;
        t_dis = 0.01;
        NP = 50;
        Max_Gen = 200;
        Rs = 0.01;
    elseif test_func_num == 6
        bounds = [0 1];
        optima = [0.079699 0.246655 0.450627 0.681420 0.933895]';
        tolerance = 0.000001;
        NP = 50;
        Max_Gen = 200;
        Rs = 0.01;
        t_dis = 0.001;
    elseif test_func_num == 7
        bounds = [0 1];
        optima = 0.079699;
        tolerance = 0.000001;
        t_dis = 0.005;
        NP = 50;
        Max_Gen = 200;
        Rs = 0.01;
    elseif test_func_num == 8
        bounds = [-6 6; -6 6];
        optima = [3 3.5844 - 3.7793 - 2.8051; 2 - 1.8481 - 3.2832 3.1313]';
        tolerance = 0.0005;
        NP = 50;
        Max_Gen = 200;
        Rs = 0.5;
        t_dis = 0.005
    elseif test_func_num == 9
        bounds = [-1.9 1.9; -1.1 1.1];
        optima = [0.089842 - 0.089842; -0.712656 0.712656]';
        tolerance = 0.000001;
        t_dis = 0.005;
        NP = 50;
        Max_Gen = 200;
        Rs = 0.5;

    elseif test_func_num == 10
        bounds = [-65.536 65.535; -65.536 65.535];

        optima = [-32, -32];
        tolerance = 0.00001;
        t_dis = 0.05;

        NP = 50;
        Max_Gen = 200;
        Rs = 0.5;
        %
    elseif test_func_num == 11
        bounds = [-10 10; -10 10];
        load optimaM22
        optima = optimaM22;

        tolerance = 0.05;
        t_dis = 0.1;
        NP = 250;
        Max_Gen = 400;
        Rs = 0.5;
    elseif test_func_num == 12
        bounds = [-10 10; -10 10; -10 10];
        load optima9
        optima = x;
        tolerance = 0.2;
        t_dis = 0.1;
        NP = 500;
        Max_Gen = 400;
        Rs = 0.5;

    elseif test_func_num == 13
        bounds = [-10 10; -10 10; -10 10; -10 10];
        load optim18
        optima = c;
        tolerance = 0.2;
        t_dis = 0.2;
        NP = 1000;
        Max_Gen = 400;
        Rs = 0.5;

    elseif test_func_num == 14
        bounds = [0.25 10];
        optima = [0.333; 0.6242; 1.1701; 2.1933; 4.1112; 7.7063];
        tolerance = 0.0001;
        t_dis = 0.1;
        NP = 100;
        Max_Gen = 200;
        Rs = 0.2;

    elseif test_func_num == 15
        bounds = [0.25 10; 0.25 10];
        load optima19
        optima = b;
        tolerance = 0.001;
        t_dis = 0.1;
        NP = 500;
        Max_Gen = 400;
        Rs = 0.2;

    elseif test_func_num == 16
        bounds = [0.25 10; 0.25 10; 0.25 10];
        load optima20
        optima = c;
        tolerance = 0.001;
        t_dis = 0.1;
        NP = 1000;
        Max_Gen = 400;

        Rs = 0.2;

    end

    if niching_method_num == 1

        for i = 1:runs,

            r2pso(test_func_num, i, bounds, optima, tolerance, NP, Max_Gen, t_dis)

        end

    elseif niching_method_num == 2

        for i = 1:runs,

            r3pso(test_func_num, i, bounds, optima, tolerance, NP, Max_Gen, t_dis)

        end

    elseif niching_method_num == 3

        for i = 1:runs,

            r2psolhc(test_func_num, i, bounds, optima, tolerance, NP, Max_Gen, t_dis)

        end

    elseif niching_method_num == 4

        for i = 1:runs,

            r3psolhc(test_func_num, i, bounds, optima, tolerance, NP, Max_Gen, t_dis)

        end

    elseif niching_method_num == 5

        for i = 1:runs,

            cde(test_func_num, i, bounds, optima, tolerance, NP, Max_Gen, t_dis)

        end

    elseif niching_method_num == 6

        for i = 1:runs,

            ncde(test_func_num, i, bounds, optima, tolerance, NP, Max_Gen, t_dis)

        end

    elseif niching_method_num == 7

        for i = 1:runs,

            sde(test_func_num, i, bounds, optima, tolerance, NP, Max_Gen, Rs, t_dis)

        end

    elseif niching_method_num == 8

        for i = 1:runs,

            ferpso(test_func_num, i, bounds, optima, tolerance, NP, Max_Gen, t_dis)

        end

    elseif niching_method_num == 9

        for i = 1:runs
            nsde(test_func_num, i, bounds, optima, tolerance, NP, Max_Gen, Rs, t_dis)
        end

    elseif niching_method_num == 10

        for i = 1:runs
            nshde(test_func_num, i, bounds, optima, tolerance, NP, Max_Gen, Rs, t_dis)
        end

    end
