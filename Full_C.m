function F2G2=Full_C(R,pos1,pos2,latticeconstant)

%Graphene Shape
a_cc = latticeconstant/sqrt(3);
a_bn2 = sqrt((latticeconstant*1.5)^2+(a_cc/2)^2);
a_bn3 = sqrt((latticeconstant*1.5)^2+(a_cc/2)^2);
a_bn4 = sqrt((latticeconstant*1.5)^2+(a_cc*3.5)^2);
a_bn5 = sqrt((latticeconstant*4)^2+(a_cc)^2);
a_bn6 = sqrt((latticeconstant*4)^2+(latticeconstant)^2);
a_bn7 = sqrt((latticeconstant*2)^2+(a_cc*5)^2);

a_g4 = sqrt((latticeconstant*2.5)^2+(a_cc*1.5)^2);
a_g7 = sqrt((latticeconstant*6)^2+(latticeconstant)^2);

posx = R(1);
posy = R(2);
pos_dis = R(4);

%Parameters
tf1 = -3.0307;
tf2 = -0.30358;
tBN3 = 0.021313;
tBN4 = -0.015464;
tBN5 = -0.023004;
tBN6 = 0.008336;
tBN7 = -0.006155;
tBN8 = 0.001795;
tBN9 = -0.001344;
tBN10 = 0.000357;

tg1 = 0.24463;
tg2 = 0.059532;
tg3 = -0.021509;
tg4 = 0.004434;
tg5 = -0.002252;
tg6 = 0.002279;
tg7 = -0.000576;

tg1b = 0.24466;
tg2b = 0.059609;
tg3b = -0.022454;
tg4b = 0.004546;
tg5b = -0.002171;
tg6b = 0.002368;
tg7b = -0.000675;

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

if pos1(5) ~= 3 || pos2(5) ~= 3
    fprintf('Atom type goes wrong in F2G2_C part\n')
end

% f1
if pos1(6) == 1 && pos2(6) == 2
    if pos_dis < a_cc+0.1*a_cc && pos_dis > a_cc-0.1*a_cc
        F2G2 = tf1*exp(-3.37*(pos_dis-acc0)/acc0);
    elseif pos_dis < a_cc*2+0.1*a_cc && pos_dis > a_cc*2-0.1*a_cc
        F2G2 = tf2*exp(-3.37*(pos_dis-acc0*2)/(acc0*2));
    elseif pos_dis < a_bn2+0.1*a_cc && pos_dis > a_bn2-0.1*a_cc
        F2G2 = tBN3*exp(-3.37*(pos_dis-r02)/r02);
    elseif pos_dis < a_bn3+0.1*a_cc && pos_dis > a_bn3-0.1*a_cc
        F2G2 = tBN4*exp(-3.37*(pos_dis-r03)/r03);
    elseif pos_dis < a_cc*4+0.1*a_cc && pos_dis > a_cc*4-0.1*a_cc
        F2G2 = tBN5*exp(-3.37*(pos_dis-acc0*4)/(acc0*4));
    elseif pos_dis < a_bn4+0.1*a_cc && pos_dis > a_bn4-0.1*a_cc
        F2G2 = tBN6*exp(-3.37*(pos_dis-r04)/r04);
    elseif pos_dis < a_cc*5+0.1*a_cc && pos_dis > a_cc*5-0.1*a_cc
        F2G2 = tBN7*exp(-3.37*(pos_dis-acc0*5)/(acc0*5));
    elseif pos_dis < a_bn5+0.1*a_cc && pos_dis > a_bn5-0.1*a_cc
        F2G2 = tBN8*exp(-3.37*(pos_dis-r05)/r05);
    elseif pos_dis < a_bn6+0.1*a_cc && pos_dis > a_bn6-0.1*a_cc
        F2G2 = tBN9*exp(-3.37*(pos_dis-r06)/r06);
    elseif pos_dis < a_bn7+0.1*a_cc && pos_dis > a_bn7-0.1*a_cc
        F2G2 = tBN10*exp(-3.37*(pos_dis-r07)/r07);
    else
        %fprintf('Did not define F2G2 model in F2G2_c\n')
        F2G2 = 0;
    end
elseif pos1(6) == 2 && pos2(6) == 1
    if pos_dis < a_cc+0.1*a_cc && pos_dis > a_cc-0.1*a_cc
        F2G2 = tf1*exp(-3.37*(pos_dis-acc0)/acc0);
    elseif pos_dis < a_cc*2+0.1*a_cc && pos_dis > a_cc*2-0.1*a_cc
        F2G2 = tf2*exp(-3.37*(pos_dis-acc0*2)/(acc0*2));
    elseif pos_dis < a_bn2+0.1*a_cc && pos_dis > a_bn2-0.1*a_cc
        F2G2 = tBN3*exp(-3.37*(pos_dis-r02)/r02);
    elseif pos_dis < a_bn3+0.1*a_cc && pos_dis > a_bn3-0.1*a_cc
        F2G2 = tBN4*exp(-3.37*(pos_dis-r03)/r03);
    elseif pos_dis < a_cc*4+0.1*a_cc && pos_dis > a_cc*4-0.1*a_cc
        F2G2 = tBN5*exp(-3.37*(pos_dis-acc0*4)/(acc0*4));
    elseif pos_dis < a_bn4+0.1*a_cc && pos_dis > a_bn4-0.1*a_cc
        F2G2 = tBN6*exp(-3.37*(pos_dis-r04)/r04);
    elseif pos_dis < a_cc*5+0.1*a_cc && pos_dis > a_cc*5-0.1*a_cc
        F2G2 = tBN7*exp(-3.37*(pos_dis-acc0*5)/(acc0*5));
    elseif pos_dis < a_bn5+0.1*a_cc && pos_dis > a_bn5-0.1*a_cc
        F2G2 = tBN8*exp(-3.37*(pos_dis-r05)/r05);
    elseif pos_dis < a_bn6+0.1*a_cc && pos_dis > a_bn6-0.1*a_cc
        F2G2 = tBN9*exp(-3.37*(pos_dis-r06)/r06);
    elseif pos_dis < a_bn7+0.1*a_cc && pos_dis > a_bn7-0.1*a_cc
        F2G2 = tBN10*exp(-3.37*(pos_dis-r07)/r07);
    else
        %fprintf('Did not define F2G2 model in F2G2_c\n')
        F2G2 = 0;
    end
elseif pos1(6) == 1 && pos2(6) == 1
    if pos_dis < latticeconstant+0.1*a_cc && pos_dis > latticeconstant-0.1*a_cc
        F2G2 = tg1*exp(-3.37*(pos_dis-r0)/r0);
    elseif pos_dis < a_cc*3+0.1*a_cc && pos_dis > a_cc*3-0.1*a_cc
        F2G2 = tg2*exp(-3.37*(pos_dis-acc0*3)/(acc0*3));
    elseif pos_dis < latticeconstant*2+0.1*a_cc && pos_dis > latticeconstant*2-0.1*a_cc
        F2G2 = tg3*exp(-3.37*(pos_dis-r0*2)/(r0*2));
    elseif pos_dis < a_g4+0.1*a_cc && pos_dis > a_g4-0.1*a_cc
        F2G2 = tg4*exp(-3.37*(pos_dis-r0g4)/r0g4);
    elseif pos_dis < latticeconstant*3+0.1*a_cc && pos_dis > latticeconstant*3-0.1*a_cc
        F2G2 = tg5*exp(-3.37*(pos_dis-r0*3)/(r0*3));
    elseif pos_dis < a_cc*6+0.1*a_cc && pos_dis > a_cc*6-0.1*a_cc
        F2G2 = tg6*exp(-3.37*(pos_dis-acc0*6)/(acc0*6));
    elseif pos_dis < a_g7+0.1*a_cc && pos_dis > a_g7-0.1*a_cc
        F2G2 = tg7*exp(-3.37*(pos_dis-r0g7)/r0g7);
    else
        %fprintf('Did not define F2G2 model in F2G2_c\n')
        F2G2 = 0;
    end
elseif pos1(6) == 2 && pos2(6) == 2
    if pos_dis < latticeconstant+0.1*a_cc && pos_dis > latticeconstant-0.1*a_cc
        F2G2 = tg1b*exp(-3.37*(pos_dis-r0)/r0);
    elseif pos_dis < a_cc*3+0.1*a_cc && pos_dis > a_cc*3-0.1*a_cc
        F2G2 = tg2b*exp(-3.37*(pos_dis-acc0*3)/(acc0*3));
    elseif pos_dis < latticeconstant*2+0.1*a_cc && pos_dis > latticeconstant*2-0.1*a_cc
        F2G2 = tg3b*exp(-3.37*(pos_dis-r0*2)/(r0*2));
    elseif pos_dis < a_g4+0.1*a_cc && pos_dis > a_g4-0.1*a_cc
        F2G2 = tg4b*exp(-3.37*(pos_dis-r0g4)/r0g4);
    elseif pos_dis < latticeconstant*3+0.1*a_cc && pos_dis > latticeconstant*3-0.1*a_cc
        F2G2 = tg5b*exp(-3.37*(pos_dis-r0*3)/(r0*3));
    elseif pos_dis < a_cc*6+0.1*a_cc && pos_dis > a_cc*6-0.1*a_cc
        F2G2 = tg6b*exp(-3.37*(pos_dis-acc0*6)/(acc0*6));
    elseif pos_dis < a_g7+0.1*a_cc && pos_dis > a_g7-0.1*a_cc
        F2G2 = tg7b*exp(-3.37*(pos_dis-r0g7)/r0g7);
    else
        %fpr8intf('Did not define F2G2 model in F2G2_c\n')
        F2G2 = 0;
    end
else
    %fprintf('Did not define F2G2 model in F2G2_c\n')
    F2G2 = 0;
end

