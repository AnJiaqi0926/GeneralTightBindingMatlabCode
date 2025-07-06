function [Omega,Chern_number,Chern_number_v1,Chern_number_v2,Ene_v] = Chern_Calculation(band,Ham_Init,disp_vec,Nk,Hsize,kx_v,ky_v,dkx,dky,b1,b2)

N_neighbor = size(disp_vec, 1);  % 近邻原胞数量
Ene_v = zeros(Nk,Nk,Hsize);
Omega = zeros(Nk,Nk);
Chern = zeros(Nk,Nk);
Chern_v1 = zeros(Nk,Nk);
Chern_v2 = zeros(Nk,Nk);

bar = waitbar(0,'读取数据中...');
temp = 0;
for index_ky = 1:Nk
    for index_kx = 1:Nk
        kx = kx_v(index_kx);
        ky = ky_v(index_ky);
        Ham = zeros(Hsize);
        Hamkx = zeros(Hsize);
        Hamky = zeros(Hsize);
        Hamkxky = zeros(Hsize);
        for k = 1:N_neighbor
            % 计算相位因子: e^{-i k·R}
            phase = exp(-1i*(kx*disp_vec(k,1) + ky*disp_vec(k,2)));
            phasekx = exp(-1i*((kx+dkx)*disp_vec(k,1) + ky*disp_vec(k,2)));
            phaseky = exp(-1i*(kx*disp_vec(k,1) + (ky+dky)*disp_vec(k,2)));
            phasekxky = exp(-1i*((kx+dkx)*disp_vec(k,1) + (ky+dky)*disp_vec(k,2)));
            
            % 添加近邻贡献
            Ham = Ham + Ham_Init(:,:,k) * phase;
            Hamkx = Hamkx + Ham_Init(:,:,k) * phasekx;
            Hamky = Hamky + Ham_Init(:,:,k) * phaseky;
            Hamkxky = Hamkxky + Ham_Init(:,:,k) * phasekxky;
        end

        Ham = (Ham+Ham')/2;
        Hamkx = (Hamkx+Hamkx')/2;
        Hamky = (Hamky+Hamky')/2;
        Hamkxky = (Hamkxky+Hamkxky')/2;
        [vecs,vals]=eig(Ham);
        % 取实部排序
        % E可以用来画3D能带
        [E,I]=sort(real(diag(vals)));
        Ene_v(index_kx,index_ky,:) = E;
        % 对矩阵的列做同样的排序
        vecs=vecs(:,I);

        [vecskx,valskx]=eig(Hamkx);
        [Ekx,Ikx]=sort(real(diag(valskx)));
        vecskx=vecskx(:,Ikx); % 略偏离kx的波函数

        [vecsky,valsky]=eig(Hamky);
        [Eky,Iky]=sort(real(diag(valsky)));
        vecsky=vecsky(:,Iky); % 略偏离ky的波函数

        [vecskxky,valskxky]=eig(Hamkxky);
        [Ekxky,Ikxky]=sort(real(diag(valskxky)));
        vecskxky=vecskxky(:,Ikxky); % 略偏离kx和ky的波函数

        % Wilson loop
        matrix1=zeros(length(band));
        matrix2=zeros(length(band));
        matrix3=zeros(length(band));
        matrix4=zeros(length(band));
        for k1=1:length(band)
            for k2=1:length(band)
                matrix1(k1,k2)=vecs(:,band(k1))'*vecskx(:,band(k2));
                matrix2(k1,k2)=vecskx(:,band(k1))'*vecskxky(:,band(k2));
                matrix3(k1,k2)=vecskxky(:,band(k1))'*vecsky(:,band(k2));
                matrix4(k1,k2)=vecsky(:,band(k1))'*vecs(:,band(k2));
            end
        end
        line1=det(matrix1);
        line2=det(matrix2);
        line3=det(matrix3);
        line4=det(matrix4);

        % 贝利曲率和陈数
        phi=-angle(line1*line2*line3*line4);
        Omega(index_kx,index_ky)=phi;
        %Chern(index_kx,index_ky)=phi/(2*pi); 
        %这里的条件判断把陈数计算限制在第一布里渊区，kx>0,kx<0的判断对应两个谷陈数
        if ky*(kx-(b1(1)+b2(1))) < kx*(ky-(b1(2)+b2(2))) && ky*(kx+b1(1)) >= kx*(ky+b1(2)) && (ky-b2(2))*(kx-(b1(1)+b2(1))) >= kx*(ky-(b1(2)+b2(2))) && (ky-b2(2))*(kx+b1(1)) < kx*(ky+b1(2))
            if kx >= 0
                Chern_v1(index_kx,index_ky)=phi/(2*pi);
            elseif kx < 0
                Chern_v2(index_kx,index_ky)=phi/(2*pi);
            end
            Chern(index_kx,index_ky) = phi/(2*pi);
        end
        temp = temp+1;
        str=['计算中...',num2str(round(100*temp/(Nk*Nk))),'%'];
        waitbar(temp/(Nk*Nk),bar,str) 
    end
end
Chern_number=sum(Chern,'all');
Chern_number_v1=sum(Chern_v1,'all');
Chern_number_v2=sum(Chern_v2,'all');
Name = 'Chern_'+string(datestr(now,'mm_dd_HHMM'))+'.mat';
save(Name,"Chern")
Name = 'Chernv1_'+string(datestr(now,'mm_dd_HHMM'))+'.mat';
save(Name,"Chern_v1")
Name = 'Chernv2_'+string(datestr(now,'mm_dd_HHMM'))+'.mat';
save(Name,"Chern_v2")
close(bar)                % 循环结束可以关闭进度条，个人一般留着不关闭
toc;                      % tic;与toc;配合使用能够返回程序运行时间
end