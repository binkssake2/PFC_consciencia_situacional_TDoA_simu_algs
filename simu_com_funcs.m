clear
%simulação
freq_clk = 168*10^6;
rms_jitter = 15*10^(-12);
dists_f = [];
areas_f = [];
passo = 50;
range1 = -1000;
range2 = 1000;
r1 = [0,0];
r2 = [-200,0];
r3 = [200,0];

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
        
        try
            angle = direction_f(X,Y);
            [a,d] = metricas_f(freq_clk,rms_jitter,X,Y,r1,r2,r3,angle);
        catch
            a = 10^5;
            d = 10^5;
        end
        
        dists_f(end+1) = d;
        areas_f(end+1) = a;
        a
        d
    end
end

[X,Y] = meshgrid(range1:passo:range2,range1:passo:range2);
D = reshape(dists_f, [size(X,2),size(X,2)]);
A = reshape(areas_f, [size(X,2),size(X,2)]);
clearvars -except A D X Y
save simu_com_funcs_reta.mat