function [a1,a2,atom_position,atom_position_Init]=Structure_GBN_ReadXYZ(Nx, Ny, plot_structure,params)
%Current version could only build 2 layer, Optimize it to including more layers
Flag = params.Flag;
layer1 = params.layer1;
layer2 = params.layer2;

if Flag == 'ReadXYZ'
    [x,y,z] = textread('generateG.xyz','%f%f%f',3);
    a1 = [x(1),y(1)];
    a2 = [x(2),y(2)];
    a3 = [x(3),y(3),z(3)];

    [index,sublattice] = textread('sublatticesSorted.dat','%f%f');
    [pos_1,pos_x,pos_y,pos_z] = textread('generateGInit.xyz','%f%f%f%f','headerlines',4);
    [pos_1c,pos_xc,pos_yc,pos_zc] = textread('generateG.xyz','%f%f%f%f','headerlines',4);
    [x,y,z,d,dx,dy] = textread('displacementsCarrFormatall.txt','%f%f%f%f%f%f');
    
    %Lattice vector
    if not(length(pos_1)==length(pos_x) && length(pos_x)==length(pos_y) && length(pos_y)==length(pos_z) && length(index)==length(pos_1)) 
        fprintf('Read Position incorrect')
        pause
    end
    if not(length(pos_1c)==length(pos_1))
        fprintf('Rigid and Relaxed structure is not correct')
        pause
    end
    atom_position_Init = [];
    atom_position = [];
    for i = 1:length(pos_1)
        if pos_z(i) == layer2
            a1i = [pos_x(i),pos_y(i),pos_z(i),1,pos_1(i),sublattice(i),dx(i),dy(i),d(i)]; % a1 = [x, y, z, layer, atomtype] # Atom type 1 = B, 2 = N, 3 = C
            a1c = [pos_xc(i),pos_yc(i),pos_zc(i),1,pos_1c(i),sublattice(i),dx(i),dy(i),d(i)];
            atom_position_Init = [atom_position_Init;a1i];
            atom_position = [atom_position;a1c];
        elseif pos_z(i) == layer1
            a1i = [pos_x(i),pos_y(i),pos_z(i),2,pos_1(i),sublattice(i),dx(i),dy(i),d(i)];
            a1c = [pos_xc(i),pos_yc(i),pos_zc(i),2,pos_1c(i),sublattice(i),dx(i),dy(i),d(i)];
            atom_position_Init = [atom_position_Init;a1i];
            atom_position = [atom_position;a1c];
        else
            fprintf('More layer in generate.xyz file, check if necessary to include.')
        end
    end
elseif Flag == 'BuildStraightly'
    a = 2.46;
    a1 = [a,0];
    a2 = [a/2,a*sqrt(3)/2];
    a = 2.5;
    a1bn = [a,0];
    a2bn = [a/2,a*sqrt(3)/2];
    pos1 = [(a1(1)+a2(1))/3,(a1(2)+a2(2))/3,19.175,1,3,1,0,0,0];
    pos2 = [2*(a1(1)+a2(1))/3,2*(a1(2)+a2(2))/3,19.175,1,3,2,0,0,0];
    pos3 = [(a1bn(1)+a2bn(1))/3,(a1bn(2)+a2bn(2))/3,15.825,2,1,3,0,0,0];
    pos4 = [2*(a1bn(1)+a2bn(1))/3,2*(a1bn(2)+a2bn(2))/3,15.825,2,2,4,0,0,0];
    atom_position = [pos1;pos2;pos3;pos4];
    atom_position_Init = atom_position;
end
