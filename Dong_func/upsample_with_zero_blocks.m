function Y = upsample_with_zero_blocks(P, scale)
    % P: input matrix
    % scale: upscaling factor (e.g., 3)
    
    [n, m] = size(P);  % assuming square, but works for rectangular too
    Y = zeros(n * scale, m * scale);
    
    for i = 1:n
        for j = 1:m
            row = (i - 1) * scale + 1;
            col = (j - 1) * scale + 1;
            Y(row, col) = P(i, j);
        end
    end
end

