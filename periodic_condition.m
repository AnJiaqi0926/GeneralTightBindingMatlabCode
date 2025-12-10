function [temp_disx,temp_disy,temp_dis]=periodic_condition(Cell,temp_disx,temp_disy,temp_dis,pos1,pos2,a1,a2)
index = 0;
for i = -Cell:Cell
    for j = -Cell:Cell
        index = index+1;
        temp_disx(index) = pos1(1)-pos2(1)+i*a1(1)+j*a2(1);        
        temp_disy(index) = pos1(2)-pos2(2)+i*a1(2)+j*a2(2);
        temp_dis(index) = sqrt(temp_disx(index)^2+temp_disy(index)^2);
    end
end
