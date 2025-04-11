clc;
clear;

Nx=100;
Ny=100;
steps=300;
a=zeros(Nx,Ny);

% Still configuration
% a(5,3)=1;
%  a(5,4)=1;
%  a(5,5)=1;

%just a test
 a(51,50)=1;
 a(52,51)=1;
 a(50,52)=1;
 a(51,52)=1;
 a(52,52)=1;

% Starting position for the glider gun
start_x = 10; % Adjust as needed
start_y = 10; % Adjust as needed

% Glider gun relative coordinates
glider_gun_coords = [
    6, 2; 6, 3; 7, 2; 7, 3;
    6, 12; 7, 12; 8, 12;
    5, 13; 9, 13;
    4, 14; 10, 14;
    4, 15; 10, 15;
    7, 16;
    5, 17; 9, 17;
    6, 18; 7, 18; 8, 18;
    7, 19;
    6, 22; 5, 22; 4, 22;
    6, 23; 5, 23; 4, 23;
    3, 24; 7, 24;
    3, 26; 2, 26; 7, 26; 8, 26;
    4, 36; 5, 36; 4, 37; 5, 37
];

% Place the glider gun on the grid
for i = 1:size(glider_gun_coords, 1)
    x = mod(start_x + glider_gun_coords(i, 1) - 1, 100) + 1; % Wrap horizontally
    y = mod(start_y + glider_gun_coords(i, 2) - 1, 100) + 1; % Wrap vertically
    grid(x, y) = 1;
end

imagesc(a)



tempa=a;
obj=VideoWriter('conway_gof','MPEG-4');
 obj.Quality=100;
 obj.FrameRate=60;
 open(obj);

for q=1:steps
    for i=1:Nx
        for j=1:Ny

            curcell=a(i,j);
            counter=count(a,i,j,Nx,Ny);

            if curcell==1
                if counter<2
                    tempa(i,j)=0;
                end
                if counter>=2 && counter<4
                    tempa(i,j)=1;
                end
                if counter>4
                    tempa(i,j)=0;
                end
            else
                if counter==3
                    tempa(i,j)=1;
                end
            end


        end
    end
a=tempa;
    imagesc(a)
     drawnow
      f=getframe(gcf);
    writeVideo(obj,f);
     
end
obj.close();


function ncount=count(a,x,y,Nx,Ny)

top=mod(y - 1 - 1, Ny)+1;
bottom=mod(y + 1 - 1, Ny)+ 1;
left=mod(x - 1 - 1, Nx) + 1;
right=mod(x + 1 - 1, Nx) + 1;

neighborhood=[a(left,y),a(left,top),a(x,top),a(right,top),a(right,y),a(right,bottom),a(x,bottom),a(left,bottom)];
ncount=sum(neighborhood);
end



