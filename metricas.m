c = physconst('LightSpeed');
freq_clk = 168*10^6; %clk do controlador

syms x y
%coor receptores (equilatero de lado l)
l = 200;
r1 = [0,0];
r2 = [l,0];
r3 = [-l,0];
recep = {r1, r2, r3};
%coor transmissor (centro de triangulo de lado 30)
trans = [-500,-350];
x_trans = trans(1);
y_trans = trans(2);

%jitter
rms_jitter = 15*10^(-12);

%multiplicador de jitter para cada BER
%10-3 6.582
%10-4 7.782
%10-5 8.834
%10-6 9.784
%10-7 10.654
%10-8 11.462
%10-9 12.218
%10-10 12.934
%10-11 13.614
%10-12 14.260
%10-13 14.882
%10-14 15.478
%10-15 16.028

multiplier = 14.260;

j_pkpk = multiplier*rms_jitter;

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

t_jitter = (abs(t_1 - t_2)*freq_clk)*j_pkpk;

delta_t = abs(t_1 - t_2) + (1/(freq_clk)) + t_jitter; %somando de resolução do clk e jitter

%1 - calc funcao para os dois ramos
f_1 = sqrt((x_recep1-x).^2+(y_recep1-y).^2) - sqrt((x_recep2-x).^2+(y_recep2-y).^2) - c*delta_t;
f1_1 = sqrt((x_recep1-x).^2+(y_recep1-y).^2) - sqrt((x_recep2-x).^2+(y_recep2-y).^2) + c*delta_t;

%%%%%%%%%segundo par de ramos da hiperbole %%%%%%%%%%%%
%2 - calc delta_t
d_1 = sqrt((x_recep1-x_trans)^2+(y_recep1-y_trans)^2);
t_1 = d_1/c;

d_2 = sqrt((x_recep3-x_trans)^2+(y_recep3-y_trans)^2);
t_2 = d_2/c;
t_jitter = (abs(t_1 - t_2)*freq_clk)*j_pkpk;

delta_t = abs(t_1 - t_2) + (1/(freq_clk)) + t_jitter;

%2 - calc funcao para os dois ramos
f_2 = sqrt((x_recep1-x).^2+(y_recep1-y).^2) - sqrt((x_recep3-x).^2+(y_recep3-y).^2) - c*delta_t;
f1_2 = sqrt((x_recep1-x).^2+(y_recep1-y).^2) - sqrt((x_recep3-x).^2+(y_recep3-y).^2) + c*delta_t;

%%%%%%%%%terceiro par de ramos da hiperbole %%%%%%%%%%%%
%3 - calc delta_t
d_1 = sqrt((x_recep3-x_trans)^2+(y_recep3-y_trans)^2);
t_1 = d_1/c;

d_2 = sqrt((x_recep2-x_trans)^2+(y_recep2-y_trans)^2);
t_2 = d_2/c;
t_jitter = (abs(t_1 - t_2)*freq_clk)*j_pkpk;

delta_t = abs(t_1 - t_2) + (1/(freq_clk)) + t_jitter;

%3 - calc funcao para os dois ramos
f_3 = sqrt((x_recep3-x).^2+(y_recep3-y).^2) - sqrt((x_recep2-x).^2+(y_recep2-y).^2) - c*delta_t;
f1_3 = sqrt((x_recep3-x).^2+(y_recep3-y).^2) - sqrt((x_recep2-x).^2+(y_recep2-y).^2) + c*delta_t;

sols = solve(f_1, y);
f1 = matlabFunction(sols(1));
f2 = matlabFunction(sols(2));

sols = solve(f_2, y);
f3 = matlabFunction(sols(1));
f4 = matlabFunction(sols(2));

sols = solve(f_3, y);
f5 = matlabFunction(sols(1));
f6 = matlabFunction(sols(2));

xx = linspace(-4000,4000,1000);

yy1 = real(f1(xx));
yy2 = real(f2(xx));
yy3 = real(f3(xx));
yy4 = real(f4(xx));
yy5 = real(f5(xx));
yy6 = real(f6(xx));

YY = {yy1,yy2,yy3,yy4,yy5,yy6};
    
inter = {};
pre_inter = {};

for k = 1:5
    for i = k+1:6
        pre_inter{end+1} = InterX([xx;cell2mat(YY(k))],[xx;cell2mat(YY(i))]);
    end
end
        
for k = 1:length(pre_inter)
    m = size(cell2mat(pre_inter(1,k)));
    if m <= 5 
        for i = 1:2:m(2)*2
            if pre_inter{1,k}(i) ~= 0 && pre_inter{1,k}(i+1) ~= 0
                    inter{end+1} = [pre_inter{1,k}(i),pre_inter{1,k}(i+1)];

            end        
        end
    end
end

%inter já são as interceções das curvas

dist = {}
for i = 1:length(inter)
    x0 = inter{i}(1);
    y0 = inter{i}(2);
    %isso é a reta do DoA, para fazer dist ponto a reta
    dist{end+1} = distance_point_line(x0,y0,tand(angle),-1,0);
end

dist_mat = cell2mat(dist);
dist_sort = sort(dist_mat);

dist_triang = [inter(find(dist_sort(1) == dist_mat)), inter(find(dist_sort(2) == dist_mat)), inter(find(dist_sort(3) == dist_mat))];

baric = [(dist_triang{1}(1)+dist_triang{2}(1)+dist_triang{3}(1))/3, (dist_triang{1}(2)+dist_triang{2}(2)+dist_triang{3}(2))/3];

%dist do baricentro ao transmissor original

d = pdist([baric(1),baric(2) ; x_trans, y_trans],'euclidean')

area = polyarea([dist_triang{1}(1),dist_triang{2}(1),dist_triang{3}(1)],[dist_triang{1}(2),dist_triang{2}(2),dist_triang{3}(2)])