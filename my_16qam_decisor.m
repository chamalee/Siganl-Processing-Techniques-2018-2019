function [Sd] = my_16qam_decisor( Sr )

%real part 
Sr_R = real(Sr);
Sr_R(Sr_R <= -2) = -3;
Sr_R(Sr_R >= 2) = 3;
Sr_R(-2 < Sr_R & Sr_R < 0) = -1;
Sr_R(0 < Sr_R & Sr_R < 2) = 1;

%imaginary part
Sr_I = imag(Sr);
Sr_I(Sr_I <= -2) = -3;
Sr_I(Sr_I >= 2) = 3;
Sr_I(-2 < Sr_I & Sr_I < 0) = -1;
Sr_I(0 <= Sr_I & Sr_I < 2) = 1;

%output
Sd = Sr_R + j*Sr_I;

end

