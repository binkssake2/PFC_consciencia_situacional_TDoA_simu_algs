clear;

dists_f = [];
areas_f = [];
passo = 50;
range1 = -1000;
range2 = 1000;

%coor recepts

for X = range1:passo:range2
    X = X+5
    
    for Y = range1:passo:range2
                
        Y = Y+5
        
        if X == 0
            X = 1;
        end
        if Y == 0
            Y = 1;
        end

        c = physconst('LightSpeed');
        freq_clk = 168*10^6; %clk do controlador
        l = 200;
        r1 = [0,0];
        r2 = [l*3,0];
        r3 = [l/2,l*sqrt(3)/2];
        recep = {r1, r2, r3};


        syms x y
        
        %coor transmissor
        x_trans = X;
        y_trans = Y;

        %jitter
        rms_jitter = 15*10^(-12);

        %calc de delta_t e delta_d dois a dois
        %pegaa coordenadas
        coord1 = cell2mat(recep(1));
        x_recep1 = coord1(1);
        y_recep1 = coord1(2);

        coord2 = cell2mat(recep(2));
        x_recep2 = coord2(1);
        y_recep2 = coord2(2);

        coord3 = cell2mat(recep(3));
        x_recep3 = coord3(1);
        y_recep3 = coord3(2);
        %%%%%%%%%primeio par de ramos da hiperbole %%%%%%%%%%%%

        %1 - calc delta_t
        d_1 = sqrt((x_recep1-x_trans)^2+(y_recep1-y_trans)^2);
        t_1 = d_1/c;

        d_2 = sqrt((x_recep2-x_trans)^2+(y_recep2-y_trans)^2);
        t_2 = d_2/c;

        n_ciclos = (abs(t_1 - t_2)*freq_clk);
        r = normrnd(0,rms_jitter,1,ceil(n_ciclos));
        t_jitter = sum(r);

        delta_t = abs(t_1 - t_2) + (1/(freq_clk)) + t_jitter; %somando de resolução do clk e jitter (ou subtraindo)

        %1 - calc funcao para os dois ramos
        f_1 = sqrt((x_recep1-x).^2+(y_recep1-y).^2) - sqrt((x_recep2-x).^2+(y_recep2-y).^2) - c*delta_t;
        %f1_1 = sqrt((x_recep1-x).^2+(y_recep1-y).^2) - sqrt((x_recep2-x).^2+(y_recep2-y).^2) + c*delta_t;

        %%%%%%%%%segundo par de ramos da hiperbole %%%%%%%%%%%%
        %2 - calc delta_t
        d_1 = sqrt((x_recep1-x_trans)^2+(y_recep1-y_trans)^2);
        t_1 = d_1/c;

        d_2 = sqrt((x_recep3-x_trans)^2+(y_recep3-y_trans)^2);
        t_2 = d_2/c;
        n_ciclos = (abs(t_1 - t_2)*freq_clk);
        r = normrnd(0,rms_jitter,1,ceil(n_ciclos));
        t_jitter = sum(r);

        delta_t = abs(t_1 - t_2) + (1/(freq_clk)) + t_jitter;

        %2 - calc funcao para os dois ramos
        f_2 = sqrt((x_recep1-x).^2+(y_recep1-y).^2) - sqrt((x_recep3-x).^2+(y_recep3-y).^2) - c*delta_t;
        %f1_2 = sqrt((x_recep1-x).^2+(y_recep1-y).^2) - sqrt((x_recep3-x).^2+(y_recep3-y).^2) + c*delta_t;

        %%%%%%%%%terceiro par de ramos da hiperbole %%%%%%%%%%%%
        %3 - calc delta_t
        d_1 = sqrt((x_recep3-x_trans)^2+(y_recep3-y_trans)^2);
        t_1 = d_1/c;

        d_2 = sqrt((x_recep2-x_trans)^2+(y_recep2-y_trans)^2);
        t_2 = d_2/c;
        n_ciclos = (abs(t_1 - t_2)*freq_clk);
        r = normrnd(0,rms_jitter,1,ceil(n_ciclos));
        t_jitter = sum(r);

        delta_t = abs(t_1 - t_2) + (1/(freq_clk)) + t_jitter;

        %3 - calc funcao para os dois ramos
        f_3 = sqrt((x_recep3-x).^2+(y_recep3-y).^2) - sqrt((x_recep2-x).^2+(y_recep2-y).^2) - c*delta_t;
        %f1_3 = sqrt((x_recep3-x).^2+(y_recep3-y).^2) - sqrt((x_recep2-x).^2+(y_recep2-y).^2) + c*delta_t;
        
        
        sols = solve(f_1, y);
        f1 = matlabFunction(sols(1));
        f2 = matlabFunction(sols(2));

        sols = solve(f_2, y);
        f3 = matlabFunction(sols(1));
        f4 = matlabFunction(sols(2));

        sols = solve(f_3, y);
        f5 = matlabFunction(sols(1));
        f6 = matlabFunction(sols(2));

        xx = linspace(-range2,range2,2*1000);

        yy1 = real(f1(xx));
        yy2 = real(f2(xx));
        yy3 = real(f3(xx));
        yy4 = real(f4(xx));
        yy5 = real(f5(xx));
        yy6 = real(f6(xx));

        YY = {yy1,yy2,yy3,yy4,yy5,yy6};

        inter1 = {};
        inter2 = {};
        inter3 = {};
        pre_inter1 = {};
        pre_inter2 = {};
        pre_inter3 = {};


        %primeiro set de inter
        pre_inter1{end+1} = InterX([xx;cell2mat(YY(1))],[xx;cell2mat(YY(3))]);
        pre_inter1{end+1} = InterX([xx;cell2mat(YY(1))],[xx;cell2mat(YY(4))]);
        pre_inter1{end+1} = InterX([xx;cell2mat(YY(2))],[xx;cell2mat(YY(3))]);
        pre_inter1{end+1} = InterX([xx;cell2mat(YY(2))],[xx;cell2mat(YY(4))]);

        %segundo set de inter
        pre_inter2{end+1} = InterX([xx;cell2mat(YY(1))],[xx;cell2mat(YY(5))]);
        pre_inter2{end+1} = InterX([xx;cell2mat(YY(1))],[xx;cell2mat(YY(6))]);
        pre_inter2{end+1} = InterX([xx;cell2mat(YY(2))],[xx;cell2mat(YY(5))]);
        pre_inter2{end+1} = InterX([xx;cell2mat(YY(2))],[xx;cell2mat(YY(6))]);

        %terceiro set de inter
        pre_inter3{end+1} = InterX([xx;cell2mat(YY(3))],[xx;cell2mat(YY(5))]);
        pre_inter3{end+1} = InterX([xx;cell2mat(YY(3))],[xx;cell2mat(YY(6))]);
        pre_inter3{end+1} = InterX([xx;cell2mat(YY(4))],[xx;cell2mat(YY(5))]);
        pre_inter3{end+1} = InterX([xx;cell2mat(YY(4))],[xx;cell2mat(YY(6))]);



        for k = 1:length(pre_inter1)
            m = size(cell2mat(pre_inter1(1,k)));
            if m <= 5 
                for i = 1:2:m(2)*2
                    if pre_inter1{1,k}(i) ~= 0 && pre_inter1{1,k}(i+1) ~= 0
                            inter1{end+1} = [pre_inter1{1,k}(i),pre_inter1{1,k}(i+1)];

                    end        
                end
            end
        end

        for k = 1:length(pre_inter2)
            m = size(cell2mat(pre_inter2(1,k)));
            if m <= 5 
                for i = 1:2:m(2)*2
                    if pre_inter2{1,k}(i) ~= 0 && pre_inter2{1,k}(i+1) ~= 0
                            inter2{end+1} = [pre_inter2{1,k}(i),pre_inter2{1,k}(i+1)];

                    end        
                end
            end
        end

        for k = 1:length(pre_inter3)
            m = size(cell2mat(pre_inter3(1,k)));
            if m <= 5 
                for i = 1:2:m(2)*2
                    if pre_inter3{1,k}(i) ~= 0 && pre_inter3{1,k}(i+1) ~= 0
                            inter3{end+1} = [pre_inter3{1,k}(i),pre_inter3{1,k}(i+1)];

                    end        
                end
            end
        end

        triangs = {};
        for k = 1:length(inter1)
            for i = 1:length(inter2)
                for j = 1:length(inter3)
                    triangs{end+1} = [inter1{k};inter2{i};inter3{j}];
                end
            end
        end

        areas = {};
        for k = 1:length(triangs)
            areas{end+1} = polyarea([triangs{k}(1),triangs{k}(2),triangs{k}(3)],[triangs{k}(4),triangs{k}(5),triangs{k}(6)]);
        end
        
        try
            
%             %By Honghao Tang DoA
%             format long %The data show that as long shaping scientific
%             doa=[atan(y_trans/x_trans)]; %Direction of arrival
%             N=200;%Snapshots
%             w=[pi/4]';%Frequency
%             M=10;%Number of array elements
%             P=length(w); %The number of signal
%             lambda=150;%Wavelength
%             d=lambda/2;%Element spacing
%             snr=20;%SNA
%             D=zeros(P,M); %To creat a matrix with P row and M column
%             for k=1:P
%             D(k,:)=exp(-j*2*pi*d*sin(doa(k))/lambda*[0:M-1]); %Assignment matrix
%             end
%             D=D';
%             xx=2*exp(j*(w*[1:N])); %Simulate signal
%             x=D*xx;
%             x=x+awgn(x,snr);%Insert Gaussian white noise
%             R=x*x'; %Data covarivance matrix
%             [N,V]=eig(R); %Find the eigenvalues and eigenvectors of R
%             NN=N(:,1:M-P); %Estimate noise subspace
%             theta=-90:0.5:90; %Peak search
%             for ii=1:length(theta)
%             SS=zeros(1,length(M));
%             for jj=0:M-1
%              SS(1+jj)=exp(-j*2*jj*pi*d*sin(theta(ii)/180*pi)/lambda);
%             end
%             PP=SS*NN*NN'*SS';
%             Pmusic(ii)=abs(1/ PP);
%             end
%             Pmusic=10*log10(Pmusic/max(Pmusic)); %Spatial spectrum function
%             %plot(theta,Pmusic,'-k')
% 
%             [M,I] = max(Pmusic);
            angle = direction_f(y_trans,x_trans)
            %end do DoA
            
            areas_sort = {};
            areas_mat = {};
            areas_mat = cell2mat(areas);
            areas_sort = sort(areas_mat);
            
            k = find(areas_sort(1) == areas_mat);
            area = areas_sort(1);
            index_areas = [];
            choices = {};
            
            len_max = length(areas_sort);
            if len_max >= 5
                len_max = 5;
            end
            
            for i = 1:len_max
                k = find(areas_sort(i) == areas_mat);
                for h = 1:length(k)
                    index_areas(end+1) = k(h);
                    choices{end+1} = triangs(k(h));
                end
            end
            %choices = {};
            %k = find(areas_sort(1) == areas_mat);
            %area = areas_sort(1) 
            %for i = 1:length(k)
            %    choices{end+1} = triangs(k(i));
            %end

            dists = {};
            for k = 1:length(choices)
                baric = [(choices{k}{1}(1)+choices{k}{1}(2)+choices{k}{1}(3))/3, (choices{k}{1}(4)+choices{k}{1}(5)+choices{k}{1}(6))/3];
                dists{end+1} = distance_point_line(baric(1),baric(2),tand(angle),-1,0);
            end

            dists_sort = {};
            dists_mat = {};
            dists_mat = cell2mat(dists);
            dists_sort = sort(dists_mat);

            k = find(dists_sort(1) == dists_mat);
            k = k(1);
            baric = [(choices{k}{1}(1)+choices{k}{1}(2)+choices{k}{1}(3))/3, (choices{k}{1}(4)+choices{k}{1}(5)+choices{k}{1}(6))/3];
            d = pdist([baric(1),baric(2) ; x_trans, y_trans],'euclidean')
            %dist = dists_sort(1)
            area = areas_mat(index_areas(k))
            
            triangulin = choices(k);
        
        catch
            d = 10^5;
            area = 10^5;
        end
        
        
        dists_f(end+1) = d;
        areas_f(end+1) = area;
        
        
        clearvars -except dists_f areas_f X Y range1 range2 passo     
        
    end
end

[X,Y] = meshgrid(range1:passo:range2,range1:passo:range2);
D = reshape(dists_f, [size(X,2),size(X,2)]);
A = reshape(areas_f, [size(X,2),size(X,2)]);

save teste_last_triangs.mat
figure
surf(X,Y,A)
figure
surf(X,Y,D)
figure
surf(X,Y,A,D)
figure
surf(X,Y,D,A)

