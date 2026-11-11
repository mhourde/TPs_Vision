clear all;
close all;

% Lecture des images
Im1 = imread('../images/cancale1.pgm');
Im2 = imread('../images/cancale2.pgm');
Im3 = imread('../images/cancale3.pgm');
Im1_coul = imread('../images/cancale1.jpg');
Im2_coul = imread('../images/cancale2.jpg');
Im3_coul = imread('../images/cancale3.jpg');
 
% Affichage des deux premières images en niveaux de gris
figure;
affichage_image(Im1,'Image 1',1,2,1);
affichage_image(Im2,'Image 2',1,2,2);

% Choix des parametres
TailleFenetre = 15;
NbPoints = 100; 
k = 0.05;
seuil = 0.85;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detection des points d'interet avec Harris %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A DECOMMENTER POUR OBSERVER LA DETECTION DE HARRIS
[XY_1,Res_1] = harris(Im1,TailleFenetre,NbPoints,k);
[XY_2,Res_2] = harris(Im2,TailleFenetre,NbPoints,k);
figure;
affichage_POI(Im1,XY_1,'POI Image 1',1,2,1);
affichage_POI(Im2,XY_2,'POI Image 2',1,2,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Appariement des points d'interet %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A DECOMMENTER POUR OBSERVER LA MISE EN CORRESPONDANCE 
% avec/sans verification des contraintes
[XY_C1,XY_C2] = apparier_POI(Im1,XY_1,Im2,XY_2,TailleFenetre,seuil);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation (et verification) de l'homographie %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A DECOMMENTER QUAND HOMOGRAPHIE AURA ETE COMPLETEE
H = homographie(XY_C1,XY_C2);

%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul de la mosaique %
%%%%%%%%%%%%%%%%%%%%%%%%%
% A DECOMMENTER QUAND HOMOGRAPHIE AURA ETE VALIDEE
Imos = mosaique(Im1,Im2,H);
figure; 
affichage_image(uint8(Imos),'Mosaique obtenue a partir des 2 images initiales',1,1,1);
% SAUVEGARDE DE LA MOSAIQUE A DEUX IMAGES
imwrite(uint8(Imos),'mosaique2.pgm');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 2 pour la reconstruction                %
% A DECOMMENTER QUAND MOSAIQUEBIS AURA ETE ECRITE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Imosbis = mosaiquebis(Im1,Im2,H);
figure; 
affichage_image(uint8(Imosbis),'Mosaique obtenue a partir des 2 images initiales (version 2)',1,1,1);
% SAUVEGARDE DE LA MOSAIQUE A DEUX IMAGES VERSION 2
imwrite(uint8(Imosbis),'mosaique2_bis.pgm');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 3 pour la reconstruction avec les couleurs R, G et B %
% A DECOMMENTER QUAND MOSAIQUECOUL AURA ETE ECRITE             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Imoscoul = mosaiquecoul(Im1_coul,Im2_coul,H);
figure;
affichage_image(uint8(Imoscoul),'Mosaique obtenue a partir des 2 images couleur initiales (version 2)',1,1,1);
% SAUVEGARDE DE LA MOSAIQUE A DEUX IMAGES EN COULEUR VERSION 2
imwrite(uint8(Imoscoul),'mosaique2_coul.pgm');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 4 pour la reconstruction avec 3 images %  
% en couleurs et/ou en niveaux de gris           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detection et appariement entre Im2 et Im3
[XY_3, Res_3] = harris(Im3, TailleFenetre, NbPoints, k);
[XY_C2_3, XY_C3] = apparier_POI(Im2, XY_2, Im3, XY_3, TailleFenetre, seuil);

% Estimation de l'homographie entre Im2 et Im3
H2 = homographie(XY_C2_3, XY_C3);

% Mosaique a 3 images en niveaux de gris
Imos3 = mosaique3(Im1, Im2, Im3, H, H2);
figure; 
affichage_image(uint8(Imos3), 'Mosaique obtenue a partir des 3 images initiales', 1, 1, 1);
imwrite(uint8(Imos3), 'mosaique3.pgm');

% Mosaique a 3 images en couleur
Imos3_coul = mosaique3(Im1_coul, Im2_coul, Im3_coul, H, H2);
figure;
affichage_image(uint8(Imos3_coul), 'Mosaique couleur obtenue a partir des 3 images initiales', 1, 1, 1);
imwrite(uint8(Imos3_coul), 'mosaique3_coul.jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 5 pour la reconstruction avec N images %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Appel de mosaiqueN avec uniquement la liste d'images en entree
% (les homographies et appariements sont geres automatiquement)
images_coul = {Im1_coul, Im2_coul, Im3_coul};
ImosN = mosaiqueN(images_coul);
figure;
affichage_image(uint8(ImosN), 'Mosaique obtenue avec mosaiqueN (N images)', 1, 1, 1);
