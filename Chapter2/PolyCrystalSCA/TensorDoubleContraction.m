function [ c ] = TensorDoubleContraction( a, b )

% TensorDoubleContraction: Generates double scalar product between matrices representing tensors in Kelvin notation.
% Priyanka Dutta, 2017

c = zeros(6,6);

for i = 1:6
    for j = 1:6
        for m = 1:6
            c(i,j) = c(i,j) + a(i,m)*b(m,j);
        end
    end
end

end

