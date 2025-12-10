function F2G2=F2G2_C(R,pos1,pos2,latticeconstant)

%Graphene Shape
a_cc = latticeconstant/sqrt(3);
a_bn2 = sqrt((latticeconstant*1.5)^2+(a_cc/2)^2);

a_g4 = sqrt((latticeconstant*2.5)^2+(a_cc*1.5)^2);
a_g7 = sqrt((latticeconstant*6)^2+(latticeconstant)^2);

posx = R(1);
posy = R(2);
pos_dis = R(4);

%Parameters
tf1 = -3.0307;
tf2 = -0.19195;

tg1 = 0.24463;
tg2 = 0.066456;

tg1b = 0.24466;
tg2b = 0.067162;

if pos1(5) ~= 3 || pos2(5) ~= 3
    fprintf('Atom type goes wrong in F2G2_C part\n')
end

% f1
if pos1(6) == 1 && pos2(6) == 2
    if pos_dis < a_cc+0.1*a_cc && pos_dis > a_cc-0.1*a_cc
        F2G2 = tf1*exp(-3.37*(pos_dis-a_cc)/a_cc);
    elseif pos_dis < a_cc*2+0.1*a_cc && pos_dis > a_cc*2-0.1*a_cc
        F2G2 = tf2*exp(-3.37*(pos_dis-a_cc*2)/(a_cc*2));
    elseif pos_dis < a_bn2+0.1*a_cc && pos_dis > a_bn2-0.1*a_cc
        F2G2 = 0; %tBN3*exp(-3.37*(pos_dis-a_bn2)/a_bn2);
    else
        %fprintf('Did not define F2G2 model in F2G2_c\n')
        F2G2 = 0;
    end
elseif pos1(6) == 2 && pos2(6) == 1
    if pos_dis < a_cc+0.1*a_cc && pos_dis > a_cc-0.1*a_cc
        F2G2 = tf1*exp(-3.37*(pos_dis-a_cc)/a_cc);
    elseif pos_dis < a_cc*2+0.1*a_cc && pos_dis > a_cc*2-0.1*a_cc
        F2G2 = tf2*exp(-3.37*(pos_dis-a_cc*2)/(a_cc*2));
    elseif pos_dis < a_bn2+0.1*a_cc && pos_dis > a_bn2-0.1*a_cc
        F2G2 = 0; %tBN3*exp(-3.37*(pos_dis-a_bn2)/a_bn2);
    else
        %fprintf('Did not define F2G2 model in F2G2_c\n')
        F2G2 = 0;
    end
elseif pos1(6) == 1 && pos2(6) == 1
    if pos_dis < latticeconstant+0.1*a_cc && pos_dis > latticeconstant-0.1*a_cc
        F2G2 = tg1*exp(-3.37*(pos_dis-latticeconstant)/latticeconstant);
    elseif pos_dis < a_cc*3+0.1*a_cc && pos_dis > a_cc*3-0.1*a_cc
        F2G2 = tg2*exp(-3.37*(pos_dis-a_cc*3)/(a_cc*3));
    else
        %fprintf('Did not define F2G2 model in F2G2_c\n')
        F2G2 = 0;
    end
elseif pos1(6) == 2 && pos2(6) == 2
    if pos_dis < latticeconstant+0.1*a_cc && pos_dis > latticeconstant-0.1*a_cc
        F2G2 = tg1b*exp(-3.37*(pos_dis-latticeconstant)/latticeconstant);
    elseif pos_dis < a_cc*3+0.1*a_cc && pos_dis > a_cc*3-0.1*a_cc
        F2G2 = tg2b*exp(-3.37*(pos_dis-a_cc*3)/(a_cc*3));
    else
        %fpr8intf('Did not define F2G2 model in F2G2_c\n')
        F2G2 = 0;
    end
else
    %fprintf('Did not define F2G2 model in F2G2_c\n')
    F2G2 = 0;
end
