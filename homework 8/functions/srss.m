function vector = srss(matrix)
%SRSS Use the SRSS algorithm to combine columns of a matrix
%     Return a column vector with the same number of rows as @matrix
    vector = sqrt(sum(matrix.^2,2));
end

