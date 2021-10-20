c = physconst('LightSpeed');
C = {'g','b','r','b'}; %cores dos graphs
c_n = 1; %aux de cores
freq_clk = 168*10^6; %clk do controlador
%coor receptores
l = 200;
r1 = [0,0];
r2 = [l,0];
r3 = [-l,0];
recep = {r1, r2, r3};
%coor transmissor
trans = [200,-600];
x_trans = trans(1);
y_trans = trans(2);
%calc de delta_t e delta_d dois a dois
N = length(recep);
%jitter RMS
rms_jitter = 15*10^(-12);

for i = 2:N
    coord = cell2mat(recep(i));
    x_recep1 = coord(1);
    y_recep1 = coord(2);
    for j = 1:N-1
        coord = cell2mat(recep(j));
        x_recep2 = coord(1);
        y_recep2 = coord(2);
        
        %calc delta_t
        d_1 = sqrt((x_recep1-x_trans)^2+(y_recep1-y_trans)^2);
        t_1 = d_1/c;
        
        d_2 = sqrt((x_recep2-x_trans)^2+(y_recep2-y_trans)^2);
        t_2 = d_2/c;
        delta_t = abs(t_1 - t_2) + (1/(freq_clk))*3;
            
        %calc funcao para os dois ramos
        f = @(x,y) sqrt((x_recep1-x).^2+(y_recep1-y).^2) - sqrt((x_recep2-x).^2+(y_recep2-y).^2) - c*delta_t;
        f1 = @(x,y) sqrt((x_recep1-x).^2+(y_recep1-y).^2) - sqrt((x_recep2-x).^2+(y_recep2-y).^2) + c*delta_t;
   
        hold on
        plot(r1(1),r1(2),'r*')
        plot(r2(1),r2(2),'r*')
        plot(r3(1),r3(2),'r*')
        plot(trans(1),trans(2),'b*')
        fimplicit(f,[-100 400 -100 400],C{c_n})
        fimplicit(f1,[-100 400 -100 400],C{c_n})
        viscircles(trans,7); %plot circ por volta do transmissor
        legend({'Receptor1','Receptor2','Receptor3','Transmissor'},'Location','northwest','Orientation','vertical')
        c_n = c_n+1;
    end
end

