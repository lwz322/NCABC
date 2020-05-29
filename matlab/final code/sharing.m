function spop = sharing(pop, val, Rs)

    [sortval, sortindex] = sort(val, 'descend');
    popsort = pop(sortindex, :);
    valsort = val(sortindex);
    clear spop;
    i = 1;
    %     if size(popsort,1)>NP
    %         popsort=popsort(1:NP,:);
    %         valsort=valsort(1:NP);
    %     end

    while i <= size(popsort, 1) - 1,
        dist = zeros(size(popsort, 1), 1);
        dist(i + 1:size(popsort, 1), :) = sum((ones(size(popsort, 1) - i, 1) * popsort(i, :) - popsort(i + 1:size(popsort, 1), :)).^2, 2) < Rs.^2;
        spop(i).pop = [popsort(i, :); popsort(dist == 1, :)];
        spop(i).val = [valsort(i); valsort(dist == 1)];
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

    for i = 1:size(spop, 2)

        if size(spop(i).pop, 1) > 1

            for j = 2:size(spop(i).pop, 1)
                dis = sqrt(sum((ones(size(spop(i).pop, 1), 1) * spop(i).pop(j, :) - spop(i).pop).^2, 2));
                penalty = sum((1 - dis / Rs));
                spop(i).newval(j, :) = spop(i).val(j) / penalty;
            end

        end

        spop(i).newval(1, :) = spop(i).val(1, :);
    end

    %       for i=1:size(xx,1)
    %           dis=sqrt(sum((ones(size(xx,1),1)*xx(i,:)-xx).^2,2));
    %           chek(1:size(xx,1),:)=dis<Rs;
    %           penalty=sum(chek.*(1-dis/Rs));
    %           newval(i,:)=xxval(i)/penalty;
    %       end
end
