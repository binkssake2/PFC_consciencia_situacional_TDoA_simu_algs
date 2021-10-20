function [xpoints,ypoints] = points_from_conic(A,B,C,D,E,F, domain,number_points, noise_power)
    % Parsing parameters
    switch nargin
        case 0
            A=1;
            B=0;
            C=1;
            D=0;
            E=0;
            F=-1;
            domain=[-1,1];
            number_points=50;
            noise_power=0;
        case 7
            number_points=50;
            noise_power=0;
        case 8
            noise_power=0;
        case 9
        otherwise
            error('Wrong number of parameters. Expected 0, 6, 7 or 8 parameters.')
    end
    
    % Generating domain vector
    x = domain(1):(domain(2)-domain(1))/round(number_points/2):domain(2)-1/round(number_points/2);
    % Expression to solve
    exp = strcat(num2str(A),'*x^2+2*',num2str(B),'*x*y+',num2str(C),'*y^2+2*',...
        num2str(D),'*x+2*',num2str(E),'*y+',num2str(F),'=0');
    yy = solve(exp,'y');
    y = eval(yy);
 
    % Point generation
    xpoints=0;
    ypoints=0;
    cont=0;
    for k=1:size(y,2)
        for m=1:size(y,1)
           if ( imag(y(m,k))==0 )
               cont=cont+1;
               xpoints(cont)=x(k)    + noise_power*randn;
               ypoints(cont)=y(m,k)  + noise_power*randn;
           end
        end
    end
    
    
end
