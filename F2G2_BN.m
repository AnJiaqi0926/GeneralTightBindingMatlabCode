function F2G2=F2G2_BN(R,pos1,pos2,latticeconstant)

a_bn = latticeconstant/sqrt(3);
a_bn2 = sqrt((latticeconstant*1.5)^2+(a_bn/2)^2);
a_bn3 = sqrt((latticeconstant*1.5)^2+(a_bn/2)^2);
a_bn4 = sqrt((latticeconstant*1.5)^2+(a_bn*3.5)^2);
a_bn5 = sqrt((latticeconstant*4)^2+(a_bn)^2);
a_bn6 = sqrt((latticeconstant*4)^2+(latticeconstant)^2);
a_bn7 = sqrt((latticeconstant*2)^2+(a_bn*5)^2);

a_nn4 = sqrt((latticeconstant*2.5)^2+(a_bn*1.5)^2);
a_nn7 = sqrt((latticeconstant*6)^2+(latticeconstant)^2);

posx = R(1);
posy = R(2);
pos_dis = R(4);

%Parameters
tBN1 = -2.8928;
tBN2 = -0.15213;

tBB1 = 0.080464;
tBB2 = 0.067819;

tNN1 = 0.24597;
tNN2 = 0.048804;

if pos1(4) == 3 || pos2(4) == 3
    fprintf('Atom type goes wrong in F2G2_BN part\n')
end

% f1
if pos1(5) == 1 && pos2(5) == 2
    if pos_dis < a_bn+0.1*a_bn && pos_dis > a_bn-0.1*a_bn
        F2G2 = tBN1;
    elseif pos_dis < a_bn*2+0.1*a_bn && pos_dis > a_bn*2-0.1*a_bn
        F2G2 = tBN2;
    elseif pos_dis < a_bn2+0.1*a_bn && pos_dis > a_bn2-0.1*a_bn
        F2G2 = 0; %tBN3;
    elseif pos_dis < a_bn3+0.1*a_bn && pos_dis > a_bn3-0.1*a_bn
        F2G2 = 0; %tBN4;
    else
        F2G2 = 0;
    end
elseif pos1(5) == 2 && pos2(5) == 1
    if pos_dis < a_bn+0.1*a_bn && pos_dis > a_bn-0.1*a_bn
        F2G2 = tBN1;
    elseif pos_dis < a_bn*2+0.1*a_bn && pos_dis > a_bn*2-0.1*a_bn
        F2G2 = tBN2;
    elseif pos_dis < a_bn2+0.1*a_bn && pos_dis > a_bn2-0.1*a_bn
        F2G2 = 0; %tBN3;
    elseif pos_dis < a_bn3+0.1*a_bn && pos_dis > a_bn3-0.1*a_bn
        F2G2 = 0; %tBN4;
    else
        F2G2 = 0;
    end
elseif pos1(5) == 1 && pos2(5) == 1
    if pos_dis < latticeconstant+0.1*a_bn && pos_dis > latticeconstant-0.1*a_bn
        F2G2 = tBB1;
    elseif pos_dis < a_bn*3+0.1*a_bn && pos_dis > a_bn*3-0.1*a_bn
        F2G2 = tBB2;
    elseif pos_dis < latticeconstant*2+0.1*a_bn && pos_dis > latticeconstant*2-0.1*a_bn
        F2G2 = 0; %tBB3;
    elseif pos_dis < a_nn4+0.1*a_bn && pos_dis > a_nn4-0.1*a_bn
        F2G2 = 0; %tBB4;
    else
        F2G2 = 0;
    end
elseif pos1(5) == 2 && pos2(5) == 2
    if pos_dis < latticeconstant+0.1*a_bn && pos_dis > latticeconstant-0.1*a_bn
        F2G2 = tNN1;
    elseif pos_dis < a_bn*3+0.1*a_bn && pos_dis > a_bn*3-0.1*a_bn
        F2G2 = tNN2;
    elseif pos_dis < latticeconstant*2+0.1*a_bn && pos_dis > latticeconstant*2-0.1*a_bn
        F2G2 = 0; %tNN3;
    elseif pos_dis < a_nn4+0.1*a_bn && pos_dis > a_nn4-0.1*a_bn
        F2G2 = 0; %tNN4;
    else
        F2G2 = 0;
    end
else
    fprintf('Did not define F2G2 model in F2G2_c\n')
    F2G2 = 0;
end


