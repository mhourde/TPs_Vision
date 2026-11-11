% Calcul d'une mosaique d'image a partir de N images stockees dans un cell array : images = {I1, I2, ..., IN}.
% Seules les images sont fournies en entree : les homographies consecutives
% sont estimees automatiquement via Harris et appariement de points d'interet.
% L'image centrale k_ref = round(N/2) est choisie comme image de reference.
% Fonctionne en niveaux de gris et en couleur.

function [Imos] = mosaiqueN(images)

N = length(images);
if N < 2
    error('Il faut au moins 2 images pour construire une mosaique.');
end

% Parametres pour la detection de Harris et la mise en correspondance
TailleFenetre = 15;
NbPoints = 100;
k = 0.05;
seuil = 0.85;

% Recuperation du nombre de canaux (1 pour niveaux de gris, 3 pour couleur)
nbCanaux = size(images{1}, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETAPE 1 : Detection des points d'interet et estimation des         %
%           homographies entre images consecutives                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detection des points d'interet avec Harris sur chaque image
XY = cell(1, N);
for i = 1:N
    I_actuelle = images{i};
    if nbCanaux == 3
        I_gris = rgb2gray(I_actuelle);
    else
        I_gris = I_actuelle;
    end
    [XY{i}, ~] = harris(I_gris, TailleFenetre, NbPoints, k);
end

% Appariement et calcul des homographies consecutives
% H_consec{i} est l'homographie permettant de passer de l'image i a l'image i+1
H_consec = cell(1, N - 1);
for i = 1:(N - 1)
    I_a = images{i};
    I_b = images{i + 1};
    if nbCanaux == 3
        I_a = rgb2gray(I_a);
        I_b = rgb2gray(I_b);
    end

    [XY_Ca, XY_Cb] = apparier_POI(I_a, XY{i}, I_b, XY{i + 1}, TailleFenetre, seuil);
    H_consec{i} = homographie(XY_Ca, XY_Cb);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETAPE 2 : Composition des homographies vers l'image centrale pivot %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% On choisit l'image centrale comme repere de reference
k_ref = round(N / 2);

% H_to_ref{i} transforme les coordonnees de l'image i vers le repere de I_{k_ref}
H_to_ref = cell(1, N);
H_to_ref{k_ref} = eye(3);

% Pour les images situees a gauche de la reference (i < k_ref)
for i = (k_ref - 1):-1:1
    H_to_ref{i} = H_to_ref{i + 1} * H_consec{i};
    H_to_ref{i} = H_to_ref{i} ./ H_to_ref{i}(3, 3);
end

% Pour les images situees a droite de la reference (i > k_ref)
for i = (k_ref + 1):N
    H_inv_step = inv(H_consec{i - 1});
    H_to_ref{i} = H_to_ref{i - 1} * H_inv_step;
    H_to_ref{i} = H_to_ref{i} ./ H_to_ref{i}(3, 3);
end

% Calcul des homographies inverses pour la retro-projection (de I_{k_ref} vers l'image i)
H_from_ref = cell(1, N);
for i = 1:N
    H_inv = inv(H_to_ref{i});
    H_from_ref{i} = H_inv ./ H_inv(3, 3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETAPE 3 : Projection des coins et calcul de la taille mosaique    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_x = [];
all_y = [];
nbl = zeros(1, N);
nbc = zeros(1, N);

for i = 1:N
    [nbl(i), nbc(i), ~] = size(images{i});
    xy_coins = [1 1; nbc(i) 1; nbc(i) nbl(i); 1 nbl(i)];
    xy_coins_proj = appliquerHomographie(H_to_ref{i}, xy_coins);

    all_x = [all_x; xy_coins_proj(:, 1)];
    all_y = [all_y; xy_coins_proj(:, 2)];
end

xmin = floor(min(all_x));
xmax = ceil(max(all_x));
ymin = floor(min(all_y));
ymax = ceil(max(all_y));

nblImos = ymax - ymin + 1;
nbcImos = xmax - xmin + 1;

if nbCanaux == 3
    Imos = zeros(nblImos, nbcImos, 3);
else
    Imos = zeros(nblImos, nbcImos);
end

% Origine de l'image de reference dans le repere de la mosaique
O_ref_x = 1 - (xmin - 1);
O_ref_y = 1 - (ymin - 1);

% Conversion des images en double pour les calculs de ponderation
images_double = cell(1, N);
for i = 1:N
    images_double{i} = double(images{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETAPE 4 : Reconstruction de la mosaique et fusion ponderee        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for x = 1:nbcImos
    for y = 1:nblImos
        % Coordonnees dans le repere de l'image de reference
        x_ref = x - O_ref_x;
        y_ref = y - O_ref_y;

        val = 0;
        poids_total = 0;

        for i = 1:N
            % Coordonnees correspondantes dans l'image i
            xy_i = appliquerHomographie(H_from_ref{i}, [x_ref y_ref]);
            xi = round(xy_i(1));
            yi = round(xy_i(2));

            % Test d'appartenance a l'image i
            if (xi >= 1 && xi <= nbc(i) && yi >= 1 && yi <= nbl(i))
                % Calcul de la distance aux bords selon la position de l'image
                if i == 1
                    % Premiere image (a gauche) : distance au bord droit
                    d = nbc(1) - xi;
                elseif i == N
                    % Derniere image (a droite) : distance au bord gauche
                    d = xi;
                else
                    % Images intermediaires : distance minimale aux bords gauche et droit
                    d = min(xi, nbc(i) - xi);
                end

                if d > 0
                    val = val + d * images_double{i}(yi, xi, :);
                    poids_total = poids_total + d;
                end
            end
        end

        if poids_total > 0
            Imos(y, x, :) = val / poids_total;
        end
    end
end

end
