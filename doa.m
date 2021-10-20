fs = 8000;
t = (0:1/fs:1).';
x1 = cos(2*pi*t*300);
sULA = phased.ULA('NumElements',10,...
    'ElementSpacing',1);
sULA.Element.FrequencyRange = [100e6 300e6];
fc = 150e6;
X = collectPlaneWave(sULA,[x1],[20 0]',fc);

X = transpose(X);

d = 0.5; %distancia entre gravacoes
N = size(X); %numero de gravacoes
array = phased.ULA('NumElements',N(1),'ElementSpacing',d);
M=min(size(X)); 
%tira nivel DC
for cont=1:M
    X(cont,:)=X(cont,:)-mean(X(cont,:));    
end
%faz sinal analitico
for cont=1:M
    X(cont,:)=hilbert(X(cont,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% PARTE DO DOA%%%%%%%%%%%%%%
M=min(size(X));
nr_amostras=max(size(X));
N = 1; %numero de sinais a serem triangulados

%tira nivel DC
for cont=1:M
    X(cont,:)=X(cont,:)-mean(X(cont,:));    
end
%faz sinal analitico
for cont=1:M
    X(cont,:)=hilbert(X(cont,:));
end
%%% calculo Rx = valor esperado de (x * x_hermitiana)
%%% (somatorio de (x * x_hermitiana))/numero_de_amostras
Rx = zeros(M,M);
for n=1:nr_amostras
     aux=X(:,n);
    Rx=Rx+aux*aux';
end
Rx=Rx/nr_amostras; 
 %autovalores são sorteados em ordem decrescente, associando os N maximos
 %autovalores para o sub-espaço dos sinais e os outros definirao o sub
 %espaço do ruido. Os auto vetores correspondentes aos auto valores positi
 %vos formam o sub espaco dos sinais. 
[V,L] = eig(Rx);
%D eh a matriz dos autovalores e agora sorteando ele em ordem descrescente
[val,pos]= sort(diag(L),'descend');
%Ordenado na ordem decrescente os auto valores
V=V(:,pos);
%Tomando os D's auto vetores mais altos para formar o sub espaço dos sinais
% e os outros irao compor o auto espaco do ruido
auto_sig=V(:,1:N);
%usa-se o auto noise já que são perpendiculares ao espaço dos sinais e dará
%picos)
auto_noise=V(:,N:end);
%Fazendo os calcs dos algs 
teta=0:0.01:180;
MVDR=zeros(length(teta),1);
MUSIC=zeros(length(teta),1);

%o vetor ateta eh o vetor de referencia dado por 
%[1 exp(-1i*2*pi*d*cos(theta)) ...  exp(-1i*2*pi*d*cos(theta)(N-1)) ]. 
%para o cod exp(-1i*(0:M-1).'*pi*2*(d/lambda)*cosd(teta(cont))) considera
%que d/lambda = 1/2 para evitar o aliasing

%rodando o algoritmo
%SS* + GG* = I para os vetores ortonormais (é o que vai no denominador do
%MUSIC para que apareca os picos
for cont=1:length(teta)
    ateta=exp(-1i*(0:M-1).'*pi*cosd(teta(cont)));
    PMVDR(cont)=1/(ateta'*inv(Rx)*ateta);
    PMUSIC(cont)=1/(ateta'*(eye(M)-auto_sig*auto_sig')'*(eye(M)-auto_sig*auto_sig')*ateta);
end
teta = -90:0.01:90;
PMVDR=PMVDR/max(abs(PMVDR));
PMUSIC=PMUSIC/max(abs(PMUSIC));
figure
plot(teta,(PMUSIC),'linewidth',3)
hold on
plot(teta,(PMVDR),'--g','linewidth',2)
ylabel('Pot normalizada')
xlabel('Azimute (em graus)')
grid
legend('MVDR','MUSIC')
[M,I] = max(PMVDR);
fprintf("angulo MVDR: %i\n", round(teta(I)))
[M,I] = max(PMUSIC);
fprintf("angulo MUSIC: %i\n", round(teta(I)))