function [ imageSup ] = imageFusion( env, velEst, map, pos, velPow)

% env = env / 255;

imageSup0 = env;
imageSup1 = env;
imageSup2 = env;

velEst = velEst + 1;

[m, n] = size(velEst);

% pos(2) = pos(2) - 0.0049 / 0.1e-3;

for i = 1 : m
    for j = 1: n
        if   ((velEst(i, j) < 126) || (velEst(i, j) > 129)) && (env(pos(2) + i, pos(1) + j) < 0.17) 
            imageSup0(pos(2) + i, pos(1) + j) = map(velEst(i, j), 1);
            imageSup1(pos(2) + i, pos(1) + j) = map(velEst(i, j), 2);
            imageSup2(pos(2) + i, pos(1) + j) = map(velEst(i, j), 3);
        end
    end
end

imageSup = cat(3, imageSup0, imageSup1, imageSup2);

end

% (velEst(i, j) > 60) && (velEst(i, j) < 180) && ...