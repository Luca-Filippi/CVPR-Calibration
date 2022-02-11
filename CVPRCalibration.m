%CVPR2020, Lab Lecture #1, starter code

% Load an image containing the checkerboard pattern
imageFileName = fullfile('images','image1.tif'); 
I = imread(imageFileName);

size(I)
class(I)

figure
imshow(I)

impixelinfo


[imagePoints, boardSize] = detectCheckerboardPoints(I);
% Visualizing the detected points
% The X direction is the long side and the Y is the short side

figure
imshow(I)
hold on
for ii=1:size(imagePoints,1)
    plot(imagePoints(ii,1),imagePoints(ii,2),'or') %o for circle and r for red
    hndtxt=text(imagePoints(ii,1),imagePoints(ii,2),num2str(ii)); 
    set(hndtxt,'fontsize',18,'color','green');
end

 
clear imageData
  
iimage=[1 4 7 10]; %indices of the images to be processed
  
for ii=1:length(iimage)
      imageFileName = fullfile('images',['image' num2str(iimage(ii)) '.tif']);
  imageData(ii).I = imread(imageFileName);
  imageData(ii).XYpixel = detectCheckerboardPoints( imageData(ii).I); %#ok
  hnd=figure; %when creating a figure, a handle is returned, that can be used for accessing properties
  imshow(imageData(ii).I,'InitialMagnification',300)
  hold on
  for jj=1:size(imageData(ii).XYpixel ,1) 
    x=imageData(ii).XYpixel(jj,1); 
    y=imageData(ii).XYpixel(jj,2); 
    %plot(x,y,'or')

    hndtxt=text(x,y,num2str(jj)); 
    set(hndtxt,'fontsize',14,'color','green');
   end
pause (1)
end

%--------------------------------------------------------------------------
% Punto 1 progetto d'esame

%----Calcolo Corrispondenze----%
squaresize =30; % [mm]
for ii=1:length(iimage)
XYpixel=imageData(ii).XYpixel;
clear Xmm Ymm
for jj=1:length(XYpixel)
[row,col]=ind2sub([12,13],jj); %linear index to row,col
Xmm=(col -1)*squaresize;
Ymm=(row -1)*squaresize;

imageData(ii).XYmm(jj,:)=[Xmm Ymm];

end
pause (1)
end


for ii=1:length(iimage) 
 %----Calcolo di H----%
    XYpixel=imageData(ii).XYpixel; 
    XYmm=imageData(ii).XYmm;
    imageData(ii).A=[];
    imageData(ii).p=[];
    for jj=1:length(XYpixel)%basterebbe 4
        imageData(ii).Xpixel=XYpixel(jj,1);
        imageData(ii).Ypixel=XYpixel(jj,2);
        imageData(ii).Xmm=XYmm(jj,1);
        imageData(ii).Ymm=XYmm(jj,2);
        imageData(ii).m=[imageData(ii).Xmm; imageData(ii).Ymm; 1];%fissiamo Zmm=0
        O=[0;0;0];
        imageData(ii).A=[imageData(ii).A; imageData(ii).m' O' -imageData(ii).Xpixel*imageData(ii).m';...
            O' imageData(ii).m' -imageData(ii).Ypixel*imageData(ii).m']; %#ok 
        imageData(ii).p=[imageData(ii).p; 0; 0]; %#ok
    end
    [U,S,V]=svd(imageData(ii).A);
    imageData(ii).h=V(:,end); 
%      imageData(ii).H=[imageData(ii).h(1) imageData(ii).h(2) imageData(ii).h(3);...
%         imageData(ii).h(4) imageData(ii).h(5) imageData(ii).h(6);...
%         imageData(ii).h(7) imageData(ii).h(8) imageData(ii).h(9)]; 
      imageData(ii).H=reshape(imageData(ii).h, [3 3])';
%----Calcolo h1,h2,h3----%
    imageData(ii).h1=imageData(ii).H(:,1);
    imageData(ii).h2=imageData(ii).H(:,2);
    imageData(ii).h3=imageData(ii).H(:,3);
%---- Calcolo vij----%
    %formula lecture 9 slide 62
    imageData(ii).v11=[imageData(ii).H(1,1)*imageData(ii).H(1,1)... %prima componente
    imageData(ii).H(1,1)*imageData(ii).H(2,1)+imageData(ii).H(2,1)*imageData(ii).H(1,1)...%seconda componente
    imageData(ii).H(2,1)*imageData(ii).H(2,1)...%terza componente
    imageData(ii).H(3,1)*imageData(ii).H(1,1)+imageData(ii).H(1,1)*imageData(ii).H(3,1)...%quarta componente
    imageData(ii).H(3,1)*imageData(ii).H(2,1)+imageData(ii).H(2,1)*imageData(ii).H(3,1)...%quinta componente
    imageData(ii).H(3,1)*imageData(ii).H(3,1)];%sesta componente


    imageData(ii).v12=[imageData(ii).H(1,1)*imageData(ii).H(1,2)... %prima componente
    imageData(ii).H(1,1)*imageData(ii).H(2,2)+imageData(ii).H(2,1)*imageData(ii).H(1,2)...%seconda componente
    imageData(ii).H(2,1)*imageData(ii).H(2,2)...%terza componente
    imageData(ii).H(3,1)*imageData(ii).H(1,2)+imageData(ii).H(1,1)*imageData(ii).H(3,2)...%quarta componente
    imageData(ii).H(3,1)*imageData(ii).H(2,2)+imageData(ii).H(2,1)*imageData(ii).H(3,2)...%quinta componente
    imageData(ii).H(3,1)*imageData(ii).H(3,2)];%sesta componente

    imageData(ii).v22=[imageData(ii).H(1,2)*imageData(ii).H(1,2)... %prima componente
    imageData(ii).H(1,2)*imageData(ii).H(2,2)+imageData(ii).H(2,2)*imageData(ii).H(1,2)...%seconda componente
    imageData(ii).H(2,2)*imageData(ii).H(2,2)...%terza componente
    imageData(ii).H(3,2)*imageData(ii).H(1,2)+imageData(ii).H(1,2)*imageData(ii).H(3,2)...%quarta componente
    imageData(ii).H(3,2)*imageData(ii).H(2,2)+imageData(ii).H(2,2)*imageData(ii).H(3,2)...%quinta componente
    imageData(ii).H(3,2)*imageData(ii).H(3,2)];%sesta componente


%----Definizione matrice V----%

       imageData(ii).V=[imageData(ii).v12; (imageData(ii).v11-imageData(ii).v22);...
           imageData(ii).v12; (imageData(ii).v11-imageData(ii).v22);...
           imageData(ii).v12; (imageData(ii).v11-imageData(ii).v22)];
      [imageData(ii).U, imageData(ii).W, imageData(ii).S]=svd(imageData(ii).V);
       
%----Calcolo vettore b----%
     imageData(ii).b=imageData(ii).V(:,6);

%----Calcolo Matrici B e L----%
%imposto direttamente B33=1 dato che K33=1
     imageData(ii).B=[imageData(ii).b(1) imageData(ii).b(2) imageData(ii).b(4);...
         imageData(ii).b(2) imageData(ii).b(3) imageData(ii).b(5);...
         imageData(ii).b(4) imageData(ii).b(5) 1];
     %imageData(ii).L=chol(imageData(ii).B); %le componenti Bij sono piccole
     % per tanto Matlab non considera B definita positiva quindi non posso
     % usare la funzione chol

%----Calcolo matrice K----%

    %imageData(ii).K=inv(imageData(ii).L'); dato che non posso calcolare L
    %con chol(B) questa istruzione non va bene
    au=sqrt((1+imageData(ii).B(1,3)^2)/abs(imageData(ii).B(1,1)));
    av=sqrt((1+imageData(ii).B(2,3)^2)/abs(imageData(ii).B(2,2)));
    u0=-imageData(ii).B(1,3)*au;
    v0=-imageData(ii).B(2,3)*av;
    imageData(ii).K=[au 0 u0; 0 av v0; 0 0 1];

%----Calcolo matrice R e vettore t----%
    imageData(ii).r1=inv(imageData(ii).K)*imageData(ii).h1;
    imageData(ii).r2=inv(imageData(ii).K)*imageData(ii).h2;
    imageData(ii).r3=cross(imageData(ii).r1,imageData(ii).r2);

    imageData(ii).R=[imageData(ii).r1 imageData(ii).r2 imageData(ii).r3];

    imageData(ii).t=inv(imageData(ii).K)*imageData(ii).h3;

end
%--------------------------------------------------------------------------
%Punto 2 del prodetto d'esame
%Consideriamo imageData(1)
%---- Calcolo matrici P e Q e il vettore q----%
Q=imageData(1).K*imageData(1).R;
q=imageData(1).K*imageData(1).t;
P=[Q q];
p1=P(1,:);%prima riga di P
p2=P(2,:);
p3=P(3,:);

error=0;%errore parte da 0, man mano aumenterà per ogni iterazione del ciclo

imshow(imageData(1).I,'InitialMagnification',300)
XYpixel=imageData(1).XYpixel;

for ii=1:size(imagePoints,1)
    u=imagePoints(ii,1);
    v=imagePoints(ii,2);
    m1=[u; v; 1];%m1 corrisponde alla notazione di m' del PDF

    c1=-inv(Q)*q;
    c=[c1; 1];

    m=c+[inv(Q)*m1; 0];

    error=error+((p1*m)/(p3*m)-u)^2+((p2*m)/(p3*m)-v)^2;
    
    hndtxt=text(((p1*m)/(p3*m)),((p2*m)/(p3*m)),num2str(ii));
    set(hndtxt,'fontsize',14,'color', 'blue');
end
pause(1);

%--------------------------------------------------------------------------
%Punto 3 del progetto d'esame
for ii=1:length(iimage) 
    figure
    imshow(imageData(ii).I,'InitialMagnification',200) 
    hold on
    
    l=100; %lato quadrato

    topleftcornerX =60;
    topleftcornerY =60; 
    
    vX=topleftcornerX+[0 0 l l]; 
    vY=topleftcornerY+[0 l l 0];


    homogeneous = [vX; vY; ones(1,length(vX))];
    proj_hom=imageData(ii).H*homogeneous;
    proj=[ proj_hom(1,:)./proj_hom(3,:); proj_hom(2,:)./proj_hom(3,:); proj_hom(3,:)./proj_hom(3,:)];
    %dato che si paral di figure 3D bisogna aggiungere a proj una terza
    %dimensione (Z), in questo caso si tal lascia Z fisso a 1 per il primo
    %rettangolo e 0.9 per il secondo rettangolo.
    proj(:,end+1)=proj(:,1);

    alpha=1;
    fill3(proj(1,:),proj(2,:),proj(3,:),'red','FaceAlpha',alpha);
    
    proj2=0.9*proj;

    fill3(proj2(1,:),proj2(2,:),proj2(3,:),'green','FaceAlpha',alpha);
    hold on;
    pause(1);
end
