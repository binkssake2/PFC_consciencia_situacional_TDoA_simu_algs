load simu_com_funcs.mat

figure
surf(X,Y,A)
title('Área de abrangência (m^2) - sem tratamento');
zlabel('Área de abrangência (m^2) - sem tratamento');
xlabel('Coordenada em x(m)');
ylabel('Coordenada em y(m)');

figure
surf(X,Y,D)
title('Distância de erro(m) - sem tratamento');
zlabel('Distância de erro(m) - sem tratamento');
xlabel('Coordenada em x(m)');
ylabel('Coordenada em y(m)');



 tf1 = isoutlier(A,1);
 tf2 = isoutlier(A,2);
 tf3 = tf1+tf2;
 tf3(tf3==2) = 1;
 
 for i = 1:length(tf3)
     for k = 1:length(tf3)
         if A(i,k) == 10^5
             A(i,k) = 0;
         end
         if tf3(i,k) == 1
             A(i,k) = 0;
         end
     end
 end
 
 for i = 1:length(A)
     for k = 1:length(A)
         if A(i,k) == 0
             A(i,k) = -max(max(A));
         end
     end
 end
 
 tf11 = isoutlier(D,1);
 tf21 = isoutlier(D,2);
 tf31 = tf11+tf21;
 tf31(tf31==2) = 1;
 
 for i = 1:length(tf31)
     for k = 1:length(tf31)
         if D(i,k) == 10^5
             D(i,k) = 0;
         end
         if tf31(i,k) == 1
             D(i,k) = 0;
         end
     end
 end
 
 for i = 1:length(D)
     for k = 1:length(D)
         if D(i,k) == 0
             D(i,k) = -max(max(D));
         end
     end
 end
 
 
figure
surf(X,Y,A)
title('Área de abrangência (m^2) - com tratamento');
zlabel('Área de abrangência (m^2) - com tratamento');
xlabel('Coordenada em x(m)');
ylabel('Coordenada em y(m)');
figure
surf(X,Y,D)
title('Distância de erro (m) - com tratamento');
zlabel('Distância de erro (m) - com tratamento');
xlabel('Coordenada em x(m)');
ylabel('Coordenada em y(m)');