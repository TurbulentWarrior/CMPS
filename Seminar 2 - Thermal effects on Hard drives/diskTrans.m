%lattice init
Nx=1000;
Ny=1000;
disk=zeros(Nx,Ny);

%monteCarlo init
msteps=200;
steps=10^5;

%Thermodynamical constants
kb=1;
J=1;
Tc=2.27*J/kb;

%variable temperature
T=0.5*Tc;
 
%generating random distribution for disk shape
for i=1:Nx
for j=1:Ny

if (i-Nx/2)^2+(j-Ny/2)^2>(Nx/2.5)^2 || (i-Nx/2)^2+(j-Ny/2)^2<(Nx/10)^2
disk(i,j)=0;
else

    if rand()>0.5
    disk(i,j)=-1;

    else
        disk(i,j)=1;

    end   
end

end
end

for ms=1:msteps 
for s=1:steps

 %choose a random spin to flip
 x=randi(Nx); 
 y=randi(Ny);

 %pass the non contributing spins
 cur=disk(x,y);
    if cur==0
    continue;
    end

%neighbouring spins
top=disk(x,mod(y - 1 - 1, Ny) + 1);
bottom=disk(x,mod(y + 1 - 1, Ny) + 1);
left=disk( mod(x - 1 - 1, Nx) + 1,y);
right=disk(mod(x + 1 - 1, Nx) + 1,y);

%neighbourhood energy
H=-J*cur*(right+left+top+bottom);

%changing energy
delta=-2*H;
ex=exp(H/(kb*T));

%flip condintion
if delta<0 || rand()<ex
disk(x,y)=-disk(x,y);
end

end
%geting the image
 imagesc(disk)
 pbaspect([1 1 1])
 colormap("spring")

%just some fancy label
 str = sprintf('T = %.2f T_c', T/Tc);
 text(0.98, 0.02, str, 'Units', 'normalized', ...
             'HorizontalAlignment', 'right', ...
             'VerticalAlignment', 'bottom', ...
             'FontSize', 12, 'Color', 'w', 'FontWeight', 'bold');
 drawnow
end

