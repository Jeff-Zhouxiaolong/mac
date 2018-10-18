function [ imf_matrix ] = bemd( input_image )

% BEMD This program calculates the Bidimensional EMD of a 2-d signal using
% the process of sifting. It is dependent on the function SIFT.

tic

% Make a 'double' copy of the image

input_image_d = cast(input_image,'double');

% Assigning the initial values of 'h' (data for sifting) and 'k' (imf index)

h_func = input_image_d;
k=1;

% The process of Sifting

while(k<4)
    [imf_temp residue_temp] = sift(h_func); 
    imf_matrix(:,:,k) = imf_temp; %#ok<AGROW>
    k = k+1;
    h_func = residue_temp;
end

% Assigning the final residue to the last IMF index

imf_matrix(:,:,k) = residue_temp;

% End of BEMD Computation

toc

end
