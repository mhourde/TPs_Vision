% Calcul d'une mosaique d'image a partir de 3 images : I1, I2 et I3
% connaissant l'homographie H1 (entre I1 et I2) et H2 (entre I2 et I3).
% On choisit l'image centrale I2 comme image de reference.
% Fonctionne en niveaux de gris et en couleur.

function [Imos] = mosaique3(I1, I2, I3, H1, H2)

% On recupere la taille des trois images
[nblI1, nbcI1, nbCanaux] = size(I1);
[nblI2, nbcI2, ~] = size(I2);
[nblI3, nbcI3, ~] = size(I3);

% Normalisation des homographies
H1 = H1 ./ H1(3,3);
H2 = H2 ./ H2(3,3);

% Calcul des homographies inverses
% H1inv permet de passer de I2 vers I1
% H2inv permet de passer de I3 vers I2
H1inv = inv(H1);
H1inv = H1inv ./ H1inv(3,3);

H2inv = inv(H2);
H2inv = H2inv ./ H2inv(3,3);

% Calcul des quatre coins dans I1 et I3
xy_coinsI1 = [1 1; nbcI1 1; nbcI1 nblI1; 1 nblI1];
xy_coinsI3 = [1 1; nbcI3 1; nbcI3 nblI3; 1 nblI3];

% Projection des coins de I1 et I3 dans le repere de I2
xy_coinsI1_R2 = appliquerHomographie(H1, xy_coinsI1);
xy_coinsI3_R2 = appliquerHomographie(H2inv, xy_coinsI3);

% Determination des dimensions de l'image mosaique
xmin = min([xy_coinsI1_R2(:,1)' xy_coinsI3_R2(:,1)' 1]);
ymin = min([xy_coinsI1_R2(:,2)' xy_coinsI3_R2(:,2)' 1]);
xmax = max([xy_coinsI1_R2(:,1)' xy_coinsI3_R2(:,1)' nbcI2]);
ymax = max([xy_coinsI1_R2(:,2)' xy_coinsI3_R2(:,2)' nblI2]);

xmin = floor(xmin);
ymin = floor(ymin);
xmax = ceil(xmax);
ymax = ceil(ymax);

% Taille de la mosaique
nblImos = ymax - ymin + 1;
nbcImos = xmax - xmin + 1;

if nbCanaux == 3
    Imos = zeros(nblImos, nbcImos, 3);
else
    Imos = zeros(nblImos, nbcImos);
end

% Origine de l'image I2 dans le repere de la mosaique
O2x_Rmos = 1 - (xmin - 1);
O2y_Rmos = 1 - (ymin - 1);

% Conversion en double pour les calculs de ponderation
I1 = double(I1);
I2 = double(I2);
I3 = double(I3);

% Parcours de chaque pixel de la mosaique
for x = 1:nbcImos
    for y = 1:nblImos
        % Coordonnees correspondantes dans le repere de I2 (reference)
        x_R2 = x - O2x_Rmos;
        y_R2 = y - O2y_Rmos;

        % Test d'appartenance a I2
        appartient_I2 = (x_R2 >= 1 && x_R2 <= nbcI2 && y_R2 >= 1 && y_R2 <= nblI2);

        % Coordonnees projetees dans I1 par H1inv
        xy_R1 = appliquerHomographie(H1inv, [x_R2 y_R2]);
        x_R1 = round(xy_R1(1));
        y_R1 = round(xy_R1(2));
        appartient_I1 = (x_R1 >= 1 && x_R1 <= nbcI1 && y_R1 >= 1 && y_R1 <= nblI1);

        % Coordonnees projetees dans I3 par H2
        xy_R3 = appliquerHomographie(H2, [x_R2 y_R2]);
        x_R3 = round(xy_R3(1));
        y_R3 = round(xy_R3(2));
        appartient_I3 = (x_R3 >= 1 && x_R3 <= nbcI3 && y_R3 >= 1 && y_R3 <= nblI3);

        % Calcul des poids fondes sur la distance aux bords
        val = 0;
        poids_total = 0;

        if appartient_I1
            d1 = nbcI1 - x_R1;
            if d1 > 0
                val = val + d1 * I1(y_R1, x_R1, :);
                poids_total = poids_total + d1;
            end
        end

        if appartient_I2
            % Pour I2 au centre : distance au bord gauche pour le raccord avec I1,
            % distance au bord droit pour le raccord avec I3
            d2 = min(x_R2, nbcI2 - x_R2);
            if d2 > 0
                val = val + d2 * I2(y_R2, x_R2, :);
                poids_total = poids_total + d2;
            end
        end

        if appartient_I3
            d3 = x_R3;
            if d3 > 0
                val = val + d3 * I3(y_R3, x_R3, :);
                poids_total = poids_total + d3;
            end
        end

        if poids_total > 0
            Imos(y, x, :) = val / poids_total;
        end
    end
end

end
