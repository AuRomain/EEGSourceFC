function [M, sets] = spectral_profile_summary(R, opt, elec)
%SPECTRAL_PROFILE_SUMMARY Spectral measures of every recording from its spectra and specparam fits (spectral_profile_features).
%   [M, sets] = spectral_profile_summary(R, opt, elec)
%   R      struct array, one per recording: labels, iaf (individual_alpha_frequency),
%          freq and pow (channels x frequencies), freq_fit and aperiodic
%          (the specparam aperiodic fit, channels x freq_fit), fooof
%          (per channel: offset, exponent, alpha_cf, alpha_pw (NaN without
%          an alpha peak), r2)
%   opt    S.iaf: channels (posterior), search (the alpha range), bins, bin_hz
%   elec   the cap (get_ftElec of S.montages{1}; optional): a cap without
%          the 10-05 names of a set gets the channels nearest to their
%          positions (cap_channels)
%   M      one row per recording, the measures of spectral_profile_features
%   sets   the channels of each set present in every recording (temporal:
%          FT7 FT8 T7 T8 TP7 TP8 FT9 FT10 T9 T10 TP9 TP10, or T3 T4), by
%          name, or by position on a cap with other names (cap_channels)
if nargin < 3, elec = []; end
common = R(1).labels(:)';
for i = 2:numel(R), common = intersect(common, R(i).labels(:)', 'stable'); end
sets = struct('post', {cap_channels(opt.channels, common, elec)}, ...
    'front', {cap_channels({'F3', 'F1', 'Fz', 'F2', 'F4'}, common, elec)}, ...
    'fp', {cap_channels({'Fp1', 'Fpz', 'Fp2'}, common, elec)}, ...
    'temp', {cap_channels({'FT7', 'FT8', 'T7', 'T8', 'TP7', 'TP8', 'FT9', 'FT10', 'T9', 'T10', 'TP9', 'TP10', 'T3', 'T4'}, ...
    common, elec)}, ...   % T3, T4: the old names of T7, T8 (kept by the 10-05 cap, which has both)
    'all', {common});
names = spectral_profile_measures();
n = numel(R);
M = array2table(nan(n, numel(names)), 'VariableNames', names);
M.IAF_source = strings(n, 1);  M.IAF_flag = strings(n, 1);
lo = opt.search(1);  hi = opt.search(2);
half = opt.bins * opt.bin_hz;                                % the IAF band: IAF +/- half
med = @(v) median(v, 'omitnan');
for i = 1:n
    r = R(i);
    in = @(set) ismember(r.labels, sets.(set));
    f = r.freq(:)';  ff = r.freq_fit(:)';
    [~, fi] = ismember(round(ff * 1e3), round(f * 1e3));    % the fit's frequencies in the spectrum
    P = double(r.pow);  A = double(r.aperiodic);  Pf = P(:, fi);
    q = r.fooof;
    M.IAF(i) = r.iaf.value;  M.IAF_PAF(i) = r.iaf.paf;  M.IAF_CoG(i) = r.iaf.cog;
    M.IAF_cov(i) = r.iaf.value;
    if isnan(M.IAF_cov(i)), M.IAF_cov(i) = r.iaf.cog; end
    M.IAF_IRASA(i) = r.iaf.irasa;  M.IAF_specparam(i) = r.iaf.specparam;
    M.IAF_source(i) = string(r.iaf.source);  M.IAF_flag(i) = string(r.iaf.flag);
    pw = q.alpha_pw(:);  pw0 = pw;  pw0(isnan(pw0)) = 0;
    M.alpha_pw_post(i) = med(pw0(in('post')));
    M.alpha_peak_share(i) = mean(~isnan(pw(in('post'))));
    M.alpha_cf_post(i) = med(q.alpha_cf(in('post')));
    M.alpha_front_post(i) = med(pw0(in('front'))) - M.alpha_pw_post(i);
    M.exponent_all(i) = med(q.exponent(in('all')));  M.offset_all(i) = med(q.offset(in('all')));
    M.exponent_post(i) = med(q.exponent(in('post')));  M.fit_r2(i) = med(q.r2(in('all')));
    band = @(a, b) f >= a & f <= b;
    pp = mean(P(in('post'), :), 1);
    M.log_alpha_post(i) = log10(mean(pp(band(lo, hi))));
    M.rel_alpha_post(i) = sum(pp(band(lo, hi))) / sum(pp(band(1, 45)));
    M.theta_alpha_post(i) = log10(sum(pp(f >= 4 & f < lo)) / sum(pp(band(lo, hi))));
    osc = max(Pf - A, 0);                                     % rhythmic power above the 1/f fit
    po = sum(osc(in('post'), :), 1);  pt = sum(Pf(in('post'), :), 1);
    fb = @(a, b) ff >= a & ff <= b;
    M.purity_fixed(i) = sum(po(fb(lo, hi))) / sum(pt(fb(lo, hi)));
    M.capture_fixed(i) = sum(po(fb(lo, hi))) / sum(po(fb(5, 15)));
    if ~isnan(M.IAF(i))
        b = fb(M.IAF(i) - half, M.IAF(i) + half);
        M.purity_iaf(i) = sum(po(b)) / sum(pt(b));
        M.capture_iaf(i) = sum(po(b)) / sum(po(fb(5, 15)));
    end
    pT = mean(P(in('temp'), :), 1);
    M.hf_temporal(i) = log10(sum(pT(band(30, 45))) / sum(pT(band(1, 45))));
    M.eye_frontopolar(i) = log10(mean(mean(P(in('fp'), band(1, 3)))) / mean(mean(P(in('post'), band(1, 3)))));
end
end
