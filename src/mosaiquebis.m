% Calcul d'une mosaique d'image à partir de 2 images : I1 et I2
% connaissant l'homographie H entre les 2.
% Version utilisant les poids fixés par les distances au bord
function [Imos] = mosaiquebis(I1,I2,H)

% On recupere la taille des deux images. 
[nblI1 nbcI1] = size(I1);
[nblI2 nbcI2] = size(I2);

% On calcule l'homographie inverse, normalisee, 
% pour nous permettre d'effectuer la transformation de I2 vers I1. 
Hinv=inv(H);
Hinv=Hinv./Hinv(3,3);

% Calcul des dimensions de la mosaique. 
% Calcul des quatre coins dans l'image 2. 
% les coins de l'image 2 sont ranges en ligne selon : 
% haut_gauche, haut_droit, bas_droit, bas gauche. 
% Lignes et colonnes sont inversees.  
xy_coinsI2_R2 = [1 1; nbcI2 1; nbcI2 nblI2; 1 nblI2];

% Application de l'homographie Hinv sur ces coins. 
% Calcul des images des coins dans I1. 
xy_coinsI2_R1 = appliquerHomographie(Hinv,xy_coinsI2_R2)

% Determination des dimensions de l'image mosaique, 
% les xmin ymin xmax ymax, ou :
%  - xmin represente la plus petite abscisse parmi les abscisses des images 
%    des coins de I2 projetes dans I1 et les coins dans I1, 
%  - etc
% Lignes et colonnes sont inversees. 

xmin=min([xy_coinsI2_R1(:,1)' 1]);
ymin=min([xy_coinsI2_R1(:,2)' 1]);
xmax=max([xy_coinsI2_R1(:,1)' nbcI1]);
ymax=max([xy_coinsI2_R1(:,2)' nblI1]);

% On arrondit de maniere a etre certain d'avoir les coordonnees initiales
% bien comprises dans l'image. 
xmin=floor(xmin);
ymin=floor(ymin);
xmax=ceil(xmax);
ymax=ceil(ymax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTRUCTION DE LA MOSAIQUE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calcul de la taille de la mosaique. 
% Lignes et colonnes sont inversees. 
nblImos=ymax-ymin+1;
nbcImos=xmax-xmin+1;

Imos=zeros(nblImos,nbcImos);

% Calcul de l'origine de l'image I1 dans le repere de la mosaique Imos. 
O1x_Rmos = 1-(xmin-1);
O1y_Rmos = 1-(ymin-1);

% On convertit en double pour les calculs de ponderation
I1 = double(I1);
I2 = double(I2);

% Parcours de chaque pixel de la mosaique
for x = 1:nbcImos
    for y = 1:nblImos
        % Coordonnees correspondantes dans le repere de I1
        x_R1 = x - O1x_Rmos;
        y_R1 = y - O1y_Rmos;

        % Test d'appartenance a I1 (coordonnees entieres exactes)
        appartient_I1 = (x_R1 >= 1 && x_R1 <= nbcI1 && y_R1 >= 1 && y_R1 <= nblI1);

        % Coordonnees projetees dans I2 par l'homographie H
        xy_R2 = appliquerHomographie(H, [x_R1 y_R1]);
        x_R2 = round(xy_R2(1));
        y_R2 = round(xy_R2(2));

        % Test d'appartenance a I2
        appartient_I2 = (x_R2 >= 1 && x_R2 <= nbcI2 && y_R2 >= 1 && y_R2 <= nblI2);

        if appartient_I1 && appartient_I2
            % Zone de recouvrement : calcul des distances aux bords
            % d1 : distance au bord droit de I1
            % d2 : distance au bord gauche de I2
            d1 = nbcI1 - x_R1;
            d2 = x_R2;

            p1 = d1 / (d1 + d2);
            p2 = d2 / (d1 + d2);

            Imos(y, x) = p1 * I1(y_R1, x_R1) + p2 * I2(y_R2, x_R2);

        elseif appartient_I1
            Imos(y, x) = I1(y_R1, x_R1);

        elseif appartient_I2
            Imos(y, x) = I2(y_R2, x_R2);
        end
    end
end
end