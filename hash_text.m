function h = hash_text(txt, no_java)
%HASH_TEXT Fingerprint of a text: 16 hexadecimal characters.
%   h = hash_text(txt). MD5 (Java) when Java is available; otherwise two
%   32-bit FNV-1a hashes (forward and backward over the bytes), computed
%   exactly in double precision. The earlier fallback multiplied in uint64,
%   which saturates, so every text got the same key (V3-12).
%   hash_text(txt, true) uses the fallback (tests).
bytes = unicode2native(char(txt), 'UTF-8');
try
    if nargin > 1 && no_java, error('hash_text:fallback', 'fallback requested'); end
    md = java.security.MessageDigest.getInstance('MD5');
    b = typecast(md.digest(bytes(:)), 'uint8');
    h = lower(reshape(dec2hex(b, 2)', 1, []));
catch
    d = double(bytes(:)');
    h = [fnv1a(d), fnv1a(fliplr(d))];
end
h = h(1:16);
end

function s = fnv1a(bytes)
% 32-bit FNV-1a; the product with the FNV prime 16777619 = 2^24 + 403 is
% split so that every intermediate value stays below 2^53 (exact in double)
h = 2166136261;
for b = bytes
    h = bitxor(h, b);
    h = mod(mod(h, 256) * 16777216 + h * 403, 4294967296);
end
s = lower(dec2hex(h, 8));
end
