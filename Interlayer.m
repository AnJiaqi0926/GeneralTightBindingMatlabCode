function Interlayer=Interlayer(R,pos1,pos2)

posx = R(1);
posy = R(2);
posz = abs(R(3));
pos_dis = R(4);

if pos1(4) == pos2(4)
    fprintf('Layer Hopping goes wrong\n')
end

if (pos1(4) == 1 && pos2(4) == 2) && (pos1(6) == 1 && pos2(6) == 3)
    Interlayer = Prefactor(posz,3.047657,1.02796)*Hopping(R);
elseif (pos1(4) == 2 && pos2(4) == 1) && (pos1(6) == 3 && pos2(6) == 1)
    Interlayer = Prefactor(posz,3.047657,1.02796)*Hopping(R);
elseif (pos1(4) == 1 && pos2(4) == 2) && (pos1(6) == 1 && pos2(6) == 4)
    Interlayer = Prefactor(posz,3.47553,1.19218)*Hopping(R);
elseif (pos1(4) == 2 && pos2(4) == 1) && (pos1(6) == 4 && pos2(6) == 1)
    Interlayer = Prefactor(posz,3.47553,1.19218)*Hopping(R);
elseif (pos1(4) == 1 && pos2(4) == 2) && (pos1(6) == 2 && pos2(6) == 3)
    Interlayer = Prefactor(posz,3.0646389,1.049223)*Hopping(R);
elseif (pos1(4) == 2 && pos2(4) == 1) && (pos1(6) == 3 && pos2(6) == 2)
    Interlayer = Prefactor(posz,3.0646389,1.049223)*Hopping(R);
elseif (pos1(4) == 1 && pos2(4) == 2) && (pos1(6) == 2 && pos2(6) == 4)
    Interlayer = Prefactor(posz,3.4799,1.16223)*Hopping(R);
elseif (pos1(4) == 2 && pos2(4) == 1) && (pos1(6) == 4 && pos2(6) == 2)
    Interlayer = Prefactor(posz,3.4799,1.16223)*Hopping(R);
else
    Interlayer = 1*Hopping(R);
    fprintf('Not define interlayer hopping\n')
    pos1
    pos2
    pause
end
