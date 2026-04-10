function [A_out, flipped] = enforce_kernel_polarity(A_in, A_anchor)
%ENFORCE_KERNEL_POLARITY Resolve global sign ambiguity of a kernel.
%   [A_out, flipped] = enforce_kernel_polarity(A_in, A_anchor)
%
%   A_in     - input kernel (double)
%   A_anchor - optional anchor kernel to match polarity to. If empty or not
%              provided, the sign of A_in is chosen so that its largest-magnitude
%              entry is positive.
%
%   A_out    - possibly sign-flipped version of A_in
%   flipped  - logical flag indicating whether a flip was applied

    A_out = A_in;
    flipped = false;

    if nargin < 2 || isempty(A_anchor)
        % Choose sign so that the largest-magnitude element is positive
        [~, idx] = max(abs(A_in(:)));
        score = A_in(idx);
    else
        % Choose sign to maximize alignment with anchor kernel
        score = sum(A_in(:) .* A_anchor(:));
    end

    if score < 0
        A_out = -A_in;
        flipped = true;
    end
end

