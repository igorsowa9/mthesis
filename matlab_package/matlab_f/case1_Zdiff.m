function [Zabs] = case1_Zdiff(element, value, s)
%CASE1_ZDIFF Function to calculate impedance for different value of element
%   Values of all elements are original and value of 'element' is changed.
%   Impedance is calculated with respect to new values, seen from middle
%   node of LCL filter.

%     load('Aalto_data0to150org.mat'); %% ONLY ORIGINAL DATA
    load('Aalto_data0to150mod.mat');
    eval([element '=value;']); % assign 'value' to 'element'
    Z1 = imp_parallel(s*PhReact_L, 1/(s*Tuned_C)) + s*Tr3_L;
    Z2 = imp_parallel(Z1, 1/(s*Cable150_C2)) + Cable150_R+s*Cable150_L;
    Z3 = imp_parallel(Z2, 1/(s*Cable150_C1)) + s*Tr2_L;
    Z4 = imp_parallel(Z3, 1/(s*Cable33_C2)) + Cable33_R+s*Cable33_L;
    Z5 = imp_parallel(Z4, 1/(s*Cable33_C1)) + s*Tr1_L + LCL_R2+s*LCL_L2;
    Z6 = imp_parallel(Z5, 1/(s*LCL_C));
    Z7 = LCL_R1+s*LCL_L1;
    Z = imp_parallel(Z6, Z7);
    Zabs = abs(Z);
    
end

