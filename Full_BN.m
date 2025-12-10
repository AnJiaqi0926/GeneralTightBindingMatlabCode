function F2G2=Full_BN(R,pos1,pos2,latticeconstant)

% BN Shape
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
tBN2 = -0.26486;
tBN3 = 0.049096;
tBN4 = -0.007505;
tBN5 = -0.018199;
tBN6 = -0.010793;
tBN7 = -0.00344;
tBN8 = 0.002338;
tBN9 = 0.000267;
tBN10 = 0.000371;

tBB1 = 0.080464;
tBB2 = 0.051535;
tBB3 = -0.035441;
tBB4 = 0.000507;
tBB5 = -0.002904;
tBB6 = 0.001261;
tBB7 = -0.000713;

tNN1 = 0.24597;
tNN2 = 0.040432;
tNN3 = -0.023648;
tNN4 = 0.001615;
tNN5 = -0.003276;
tNN6 = 0.000827;
tNN7 = -0.000612;

r0 = 2.438977765130281; 
acc0 = r0/sqrt(3);
r02 = sqrt((r0*1.5)^2+(acc0/2)^2);
r03 = sqrt((r0*1.5)^2+(acc0/2)^2); 
r04 = sqrt((r0*1.5)^2+(acc0*3.5)^2);
r05 = sqrt((r0*4)^2+(acc0)^2);  
r06 = sqrt((r0*4)^2+(r0)^2); 
r07 = sqrt((r0*2)^2+(acc0*5)^2);   

r0g4 = sqrt((r0*2.5)^2+(acc0*1.5)^2);
r0g7 = sqrt((r0*6)^2+(r0)^2);

if pos1(4) == 3 || pos2(4) == 3
    fprintf('Atom type goes wrong in F2G2_BN part\n')
end

% f1
if pos1(5) == 1 && pos2(5) == 2
    if pos_dis < a_bn+0.1*a_bn && pos_dis > a_bn-0.1*a_bn
        F2G2 = tBN1*exp(-3.37*(pos_dis-acc0)/acc0);
    elseif pos_dis < a_bn*2+0.1*a_bn && pos_dis > a_bn*2-0.1*a_bn
        F2G2 = tBN2*exp(-3.37*(pos_dis-acc0*2)/(acc0*2));
    elseif pos_dis < a_bn2+0.1*a_bn && pos_dis > a_bn2-0.1*a_bn
        F2G2 = tBN3*exp(-3.37*(pos_dis-r02)/r02);
    elseif pos_dis < a_bn3+0.1*a_bn && pos_dis > a_bn3-0.1*a_bn
        F2G2 = tBN4*exp(-3.37*(pos_dis-r03)/r03);
    elseif pos_dis < a_bn*4+0.1*a_bn && pos_dis > a_bn*4-0.1*a_bn
        F2G2 = tBN5*exp(-3.37*(pos_dis-acc0*4)/(acc0*4));
    elseif pos_dis < a_bn4+0.1*a_bn && pos_dis > a_bn4-0.1*a_bn
        F2G2 = tBN6*exp(-3.37*(pos_dis-r04)/r04);
    elseif pos_dis < a_bn*5+0.1*a_bn && pos_dis > a_bn*5-0.1*a_bn
        F2G2 = tBN7*exp(-3.37*(pos_dis-acc0*5)/(acc0*5));
    elseif pos_dis < a_bn5+0.1*a_bn && pos_dis > a_bn5-0.1*a_bn
        F2G2 = tBN8*exp(-3.37*(pos_dis-r05)/r05);
    elseif pos_dis < a_bn6+0.1*a_bn && pos_dis > a_bn6-0.1*a_bn
        F2G2 = tBN9*exp(-3.37*(pos_dis-r06)/r06);
    elseif pos_dis < a_bn7+0.1*a_bn && pos_dis > a_bn7-0.1*a_bn
        F2G2 = tBN10*exp(-3.37*(pos_dis-r07)/r07);
    else
        F2G2 = 0;
    end
elseif pos1(5) == 2 && pos2(5) == 1
    if pos_dis < a_bn+0.1*a_bn && pos_dis > a_bn-0.1*a_bn
        F2G2 = tBN1*exp(-3.37*(pos_dis-acc0)/acc0);
    elseif pos_dis < a_bn*2+0.1*a_bn && pos_dis > a_bn*2-0.1*a_bn
        F2G2 = tBN2*exp(-3.37*(pos_dis-acc0*2)/(acc0*2));
    elseif pos_dis < a_bn2+0.1*a_bn && pos_dis > a_bn2-0.1*a_bn
        F2G2 = tBN3*exp(-3.37*(pos_dis-r02)/r02);
    elseif pos_dis < a_bn3+0.1*a_bn && pos_dis > a_bn3-0.1*a_bn
        F2G2 = tBN4*exp(-3.37*(pos_dis-r03)/r03);
    elseif pos_dis < a_bn*4+0.1*a_bn && pos_dis > a_bn*4-0.1*a_bn
        F2G2 = tBN5*exp(-3.37*(pos_dis-acc0*4)/(acc0*4));
    elseif pos_dis < a_bn4+0.1*a_bn && pos_dis > a_bn4-0.1*a_bn
        F2G2 = tBN6*exp(-3.37*(pos_dis-r04)/r04);
    elseif pos_dis < a_bn*5+0.1*a_bn && pos_dis > a_bn*5-0.1*a_bn
        F2G2 = tBN7*exp(-3.37*(pos_dis-acc0*5)/(acc0*5));
    elseif pos_dis < a_bn5+0.1*a_bn && pos_dis > a_bn5-0.1*a_bn
        F2G2 = tBN8*exp(-3.37*(pos_dis-r05)/r05);
    elseif pos_dis < a_bn6+0.1*a_bn && pos_dis > a_bn6-0.1*a_bn
        F2G2 = tBN9*exp(-3.37*(pos_dis-r06)/r06);
    elseif pos_dis < a_bn7+0.1*a_bn && pos_dis > a_bn7-0.1*a_bn
        F2G2 = tBN10*exp(-3.37*(pos_dis-r07)/r07);
    else
        F2G2 = 0;
    end
elseif pos1(5) == 1 && pos2(5) == 1
    if pos_dis < latticeconstant+0.1*a_bn && pos_dis > latticeconstant-0.1*a_bn
        F2G2 = tBB1*exp(-3.37*(pos_dis-r0)/r0);
    elseif pos_dis < a_bn*3+0.1*a_bn && pos_dis > a_bn*3-0.1*a_bn
        F2G2 = tBB2*exp(-3.37*(pos_dis-acc0*3)/(acc0*3));
    elseif pos_dis < latticeconstant*2+0.1*a_bn && pos_dis > latticeconstant*2-0.1*a_bn
        F2G2 = tBB3*exp(-3.37*(pos_dis-r0*2)/(r0*2));
    elseif pos_dis < a_nn4+0.1*a_bn && pos_dis > a_nn4-0.1*a_bn
        F2G2 = tBB4*exp(-3.37*(pos_dis-r0g4)/r0g4);
    elseif pos_dis < latticeconstant*3+0.1*a_bn && pos_dis > latticeconstant*3-0.1*a_bn
        F2G2 = tBB5*exp(-3.37*(pos_dis-r0*3)/(r0*3));
    elseif pos_dis < a_bn*6+0.1*a_bn && pos_dis > a_bn*6-0.1*a_bn
        F2G2 = tBB6*exp(-3.37*(pos_dis-acc0*6)/(acc0*6));
    elseif pos_dis < a_nn7+0.1*a_bn && pos_dis > a_nn7-0.1*a_bn
        F2G2 = tBB7*exp(-3.37*(pos_dis-r0g7)/r0g7);
    else
        F2G2 = 0;
    end
elseif pos1(5) == 2 && pos2(5) == 2
    if pos_dis < latticeconstant+0.1*a_bn && pos_dis > latticeconstant-0.1*a_bn
        F2G2 = tNN1*exp(-3.37*(pos_dis-r0)/r0);
    elseif pos_dis < a_bn*3+0.1*a_bn && pos_dis > a_bn*3-0.1*a_bn
        F2G2 = tNN2*exp(-3.37*(pos_dis-acc0*3)/(acc0*3));
    elseif pos_dis < latticeconstant*2+0.1*a_bn && pos_dis > latticeconstant*2-0.1*a_bn
        F2G2 = tNN3*exp(-3.37*(pos_dis-r0*2)/(r0*2));
    elseif pos_dis < a_nn4+0.1*a_bn && pos_dis > a_nn4-0.1*a_bn
        F2G2 = tNN4*exp(-3.37*(pos_dis-r0g4)/r0g4);
    elseif pos_dis < latticeconstant*3+0.1*a_bn && pos_dis > latticeconstant*3-0.1*a_bn
        F2G2 = tNN5*exp(-3.37*(pos_dis-r0*3)/(r0*3));
    elseif pos_dis < a_bn*6+0.1*a_bn && pos_dis > a_bn*6-0.1*a_bn
        F2G2 = tNN6*exp(-3.37*(pos_dis-acc0*6)/(acc0*6));
    elseif pos_dis < a_nn7+0.1*a_bn && pos_dis > a_nn7-0.1*a_bn
        F2G2 = tNN7*exp(-3.37*(pos_dis-r0g7)/r0g7);
    else
        F2G2 = 0;
    end
else
    fprintf('Did not define F2G2 model in F2G2_c\n')
    F2G2 = 0;
end

