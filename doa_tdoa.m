c = physconst('LightSpeed');
C = {'g','b','r','b'}; %cores dos graphs
c_n = 1; %aux de cores
freq_clk = 168*10^6; %clk do controlador

%coor receptores (equilatero de lado l)
l = 200;
r1 = [0,0];
r2 = [-l,l];
r3 = [l,0];
recep = {r1, r2, r3};
%coor transmissor (centro de triangulo de lado 30)
trans = [300,-100];
x_trans = trans(1);
y_trans = trans(2);

%jitter
rms_jitter = 15*10^(-12);


%calc de delta_t e delta_d dois a dois
N = length(recep);
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
        
        %calc n_ciclos para fazer contribuição do Jitter
        n_ciclos = (abs(t_1 - t_2)*freq_clk);
        
        %calcula contribuição do jitter, media 0, de desvios padrao = RMS
        %para cada ciclo e soma para cada um
        r = normrnd(0,rms_jitter,1,ceil(n_ciclos));
        t_jitter = sum(r);
        
        delta_t = abs(t_1 - t_2) + (1/(freq_clk)) + t_jitter; %somando resolução do clk e jitter random
        %delta_t = abs(t_1 - t_2)    
        %calc funcao para os dois ramos
        f = @(x,y) abs(sqrt((x_recep1-x).^2+(y_recep1-y).^2)) - abs(sqrt((x_recep2-x).^2+(y_recep2-y).^2)) - c*delta_t;
        f1 = @(x,y) abs(sqrt((x_recep1-x).^2+(y_recep1-y).^2)) - abs(sqrt((x_recep2-x).^2+(y_recep2-y).^2)) + c*delta_t;
           
        %funcao de plot
        hold on
        plot(r1(1),r1(2),'r*')
        plot(r2(1),r2(2),'r*')
        plot(r3(1),r3(2),'r*')
        plot(trans(1),trans(2),'b*')
        fimplicit(f,[-4000 4000 -4000  4000],C{c_n})
        fimplicit(f1,[-4000 4000 -4000 4000],C{c_n})
        %viscircles(trans,7); %plot circ por volta do transmissor
        legend({'Receptor1','Receptor2','Receptor3','Transmissor'},'Location','northwest','Orientation','vertical')
        c_n = c_n+1;
    end
end
clearvars -except  y_trans x_trans
%calc DoA em relação ao receptor em (0,0)
%By Honghao Tang
format long %The data show that as long shaping scientific
doa=[atan(y_trans/x_trans)]; %Direction of arrival
N=200;%Snapshots
w=[pi/4]';%Frequency
M=10;%Number of array elements
P=length(w); %The number of signal
lambda=150;%Wavelength
d=lambda/2;%Element spacing
snr=20;%SNA
D=zeros(P,M); %To creat a matrix with P row and M column
for k=1:P
D(k,:)=exp(-j*2*pi*d*sin(doa(k))/lambda*[0:M-1]); %Assignment matrix
end
D=D';
xx=2*exp(j*(w*[1:N])); %Simulate signal
x=D*xx;
x=x+awgn(x,snr);%Insert Gaussian white noise
R=x*x'; %Data covarivance matrix
[N,V]=eig(R); %Find the eigenvalues and eigenvectors of R
NN=N(:,1:M-P); %Estimate noise subspace
theta=-90:0.5:90; %Peak search
for ii=1:length(theta)
SS=zeros(1,length(M));
for jj=0:M-1
 SS(1+jj)=exp(-j*2*jj*pi*d*sin(theta(ii)/180*pi)/lambda);
end
PP=SS*NN*NN'*SS';
Pmusic(ii)=abs(1/ PP);
end
Pmusic=10*log10(Pmusic/max(Pmusic)); %Spatial spectrum function
[M,I] = max(Pmusic);
angle = theta(I);

if angle > -90 & angle < 90
    y = linspace(0,1000);
else 
    y = linspace(-1000,0);
end
plot(y, y*tand(angle),'--k','DisplayName','DoA - MUSIC')