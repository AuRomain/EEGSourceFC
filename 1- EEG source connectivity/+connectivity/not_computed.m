function varargout = not_computed(varargin)
%NOT_COMPUTED Connectivity arrays of a measure that could not be computed: NaN edges, zero diagonal.
%   [A, B, ...] = connectivity.not_computed(A, B, ...) with each array
%   nSignals x nSignals or nWindows x nSignals x nSignals (the windows
%   first): every off-diagonal value NaN, the diagonal 0, class and size kept.
varargout = varargin;
for k = 1:nargin
    A = varargin{k};
    n = size(A, ndims(A));
    if ismatrix(A) && size(A, 1) == n
        A(~eye(n)) = NaN;
    else
        A(:, ~eye(n)) = NaN;
    end
    varargout{k} = A;
end
end
