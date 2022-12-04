%ejercicio 5 hoja 6
f=@(x,y) 0*x;

beta=0.25;
gamma=0.5;
k=0;
c2=0.1; %el bicho de c^2

t0=0;
T=100;
dt=0.01;

Nt=T/dt;
xi = p(1,:); % primera fila de p
yi = p(2,:); % segunda fila de p
elem = t(1:3,:)';

[R M]=assema(p,t,1,1,0);
   
    A=M+c2*dt^2*R*beta;

    % necesito saber que nodos dirichlet estan en la frontera dirichlet
fron_d = find(xi == 1 | xi == -1 | yi == 1 | yi == -1 );

% quitar las filas y columnas de nodos dirichlet
A0 = A;
A0(fron_d,:) = 0;
A0(:,fron_d) = 0;
for i = fron_d
    A0(i,i) = 1;
end
% fin de quitar
   

    %uhn=sin(pi*xi').^100;
    %plot(xi,uhn)
    %title(0)
    %pause
    fi=f(xi)'; %esto no hace alta porque f(x)=0
%     uhn=0*xi';
%     vhn=xi';

uhn=0*xi';
vhn=1*(((xi-0.5).^2+(yi-0.5).^2)<0.01);
vhn=vhn';%lo queremos vertical

    vect_b=-c2*R*uhn+M*fi;
    ahn=M\vect_b;
    for n=1:Nt 
        vect_b=M*(uhn+dt*vhn+(0.5-beta)*dt^2*ahn+beta*dt^2*fi);
        vect_b(fron_d)=0;
        uu=A0\vect_b;
        %uu=[0 uu' 0]'; %vect columna

        %dudas con los uhn y ahn
%         vect_b1=1/(beta*dt^2)*M*(uhn-uhn-dt*vhn-0.5*(1-2*beta)*dt^2*ahn);
%         ahn=M\vect_b1;
        %en este caso para ahn y vhn no hace alta invertir matriz
        aa=(uu-uhn-dt*vhn-(.5-beta)*dt^2*ahn)/(beta*dt^2);
        vv=vhn+dt*(gamma*aa+(1-gamma)*ahn);

        uhn=uu;
        vhn=vv;
        ahn=aa;

        trisurf(elem, xi, yi, uhn)
        axis([-1 1 -1 1 -1 1])
        %pause(0.1) %así se queda parado 0.1 segundos
        title(n*dt) %para ver el instante de tiempo en el que estoy
        view(2)
        shading interp
        caxis([-0.1 .1])
        colorbar
        pause

    end
    

