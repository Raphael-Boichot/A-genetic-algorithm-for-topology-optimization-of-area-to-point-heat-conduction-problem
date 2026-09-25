% https://doi.org/10.1016/j.ijthermalsci.2016.05.015
% https://github.com/Raphael-Boichot/A-genetic-algorithm-for-topology-optimization-of-area-to-point-heat-conduction-problem
%
% DROP-IN REPLACEMENT for correction_mutation.m -- SAME call signature as the original:
%
%   child2 = correction_mutation(child, kp_k0, k0, number_conductive_cells, p_mutation)
%
% Self-contained: this one file has everything it needs as local functions below. Nothing else to
% copy, nothing else on the path required. No toolbox needed (base MATLAB / Octave only).
%
% WHAT IT DOES DIFFERENTLY FROM THE ORIGINAL
%   The original repairs the cell count by adding/removing single random cells next to existing
%   material. This version instead:
%     1. finds the morphological skeleton (centreline) of the current conductive phase,
%     2. makes a few SHORT, LOCAL edits driven by that skeleton: extend an existing tip outward,
%        sprout a new short branch from a random point on the skeleton, or prune a short spur back
%        to its nearest junction,
%     3. only then does a small residual snap to hit the exact cell count (a global step, but now
%        a minor cleanup rather than the main mechanism).
%   p_mutation keeps the same meaning (roughly how many cells get touched), calibrated from
%   measured average edit sizes so additions and removals roughly balance.
%
% WHY: an earlier version of this replacement worked by smoothing a signed-distance field (blend +
% threshold), which is mathematically a curvature-flow / MBO step -- it eats thin filaments and
% fuses nearby branches (anastomosis) by construction. This version never smooths the existing
% material; every edit is a short, roughly 1-D change, and everything else in the mask is left
% exactly as it was.
%
% STATUS / CAVEATS (please read before trusting this in a long run):
%   - Invariants are verified: exact conductive-cell count, border ring untouched, only k0/kc
%     values, on hundreds of synthetic and real-population test children.
%   - It has been run inside a full selection loop (population + fitness + elitism) on a small
%     synthetic proxy (25x50 grid, a stand-in solver, ~150 generations, a few seeds). In that test
%     it did NOT clearly outperform the failed smoothing-based version, and a "true original
%     algorithm" baseline run through the same proxy collapsed to a similarly sparse shape too --
%     which points at the test harness (small grid, proxy solver, generation count/mutation
%     schedule mismatched to your real run) rather than settling the question either way.
%   - It has NOT been validated on your actual 50x100 (or larger) grid, your actual solver, or over
%     the generation counts your real runs use. Please treat this as a candidate to test against
%     your original files with matched seeds, not as a proven fix.
%
% OPTIONAL 6th ARGUMENT (all fields optional, same defaults as tested above):
%   opts.n_extend, opts.n_bud, opts.n_prune   how many of each edit to attempt per call
%   opts.ext_len, opts.bud_len   [min max] length in cells for extensions / buds (default [2 5],[2 4])
%   opts.prune_max_len   a leaf is only pruned if it reaches a branch point within this many cells
%                        of its endpoint (default 6) -- keeps this a spur trim, not a real-branch cut
%   opts.angle_jitter    random angular jitter added to the estimated direction, radians (default 0.6)
%   opts.bud_radius_frac new bud's radius as a fraction of the local radius at its origin (default 0.6)
%   opts.max_dist, opts.metric   forwarded to the internal distance transform (default 12, 'octagon')

function [child2, info] = correction_mutation(child, kp_k0, k0, number_conductive_cells, p_mutation, opts)

if nargin < 6 || isempty(opts); opts = struct(); end
N = number_conductive_cells;
ext_len    = get_opt(opts,'ext_len',[2 5]);
bud_len    = get_opt(opts,'bud_len',[2 4]);
prune_max  = get_opt(opts,'prune_max_len',6);
jitter_amp = get_opt(opts,'angle_jitter',0.6);
bud_frac   = get_opt(opts,'bud_radius_frac',0.6);
max_dist   = get_opt(opts,'max_dist',12);

[height,width,~] = size(child);
kc = k0*kp_k0;
interior = false(height,width); interior(2:height-1,2:width-1) = true;
Cmask = interior & (child==kc);
Nmask = interior & (child==k0);
nC = nnz(Cmask);

M = max(1, round(p_mutation*max(N,1)));           % same "budget" definition as the original
% empirical average cell-count delta of a single extend / bud / prune (measured on 50x100
% topologies, ext_len=[2 5], bud_len=[2 4], prune_max_len=6): used only to pick, by default, how
% many of each edit to attempt so that (a) additions and removals roughly BALANCE (net change near
% zero, so the final snap-to-N step below only mops up a small residual instead of doing the bulk
% of the work) and (b) the total gross change scales with the mutation budget M. Override with
% opts.n_extend/n_bud/n_prune directly if your grid/length settings differ a lot from these.
avg_extend = 17; avg_bud = 4; avg_prune = 11;
e_default = max(1, round(M/(2*(avg_extend+avg_bud))));
p_default = max(1, round(M/(2*avg_prune)));
n_extend = get_opt(opts,'n_extend', e_default);
n_bud    = get_opt(opts,'n_bud',    e_default);
n_prune  = get_opt(opts,'n_prune',  p_default);

added = false(height,width); removed = false(height,width);
mask = Cmask;

if nC > 0
    S = skeletonize(mask);                            % the ONE full skeletonize() call for this whole
    D = local_thickness(mask);                         % operator; extend/bud update S incrementally below

    %% 1) tip extension: grow a few existing tips outward
    [ep,~,~] = skeleton_endpoints(S);
    ep = ep(randperm(numel(ep), min(n_extend, numel(ep))));
    for k = 1:numel(ep)
        path = walk_skeleton(S, ep(k), 6);
        [dr,dc] = local_direction(path, height, width);
        if isnan(dr); continue; end
        L = randi(ext_len);
        ang = atan2(dr,dc) + jitter_amp*(rand-0.5)*2;
        r0v = D(ep(k)); if r0v<=0; r0v = 1; end
        [mask, added, S] = draw_segment(mask, added, S, ep(k), ang, L, r0v, height, width, interior);
    end

    %% 2) branch budding: sprout new short branches from random points on the (now updated) skeleton
    idxS = find(S);
    if ~isempty(idxS)
        picks = idxS(randi(numel(idxS), 1, min(n_bud,numel(idxS))));
        for k = 1:numel(picks)
            [dr,dc] = local_tangent(S, picks(k), height, width);
            L = randi(bud_len);
            base_ang = atan2(dr,dc) + pi/2;           % roughly perpendicular to the local fibre
            ang = base_ang + jitter_amp*(rand-0.5)*2*pi/2;
            if rand<0.5; ang = ang + pi; end           % bud on either side at random
            r0v = max(1, round(bud_frac*max(D(picks(k)),1)));
            [mask, added, S] = draw_segment(mask, added, S, picks(k), ang, L, r0v, height, width, interior);
        end
    end

    %% 3) prune a few short spurs (on the same updated skeleton -- no re-skeletonize needed)
    [ep,~,~] = skeleton_endpoints(S);
    ep = ep(randperm(numel(ep), min(n_prune, numel(ep))));
    for k = 1:numel(ep)
        path = walk_skeleton(S, ep(k), prune_max);
        if numel(path) >= prune_max; continue; end     % didn't reach a branch point in time: not a short spur
        for j = 1:numel(path)
            r0v = D(path(j)); if r0v<=0; r0v = 1; end
            disk = octagon_disk_local(round(r0v));
            [mask, removed] = stamp(mask, removed, path(j), disk, height, width, interior, false);
            S(path(j)) = false;
        end
    end
end

child2 = child;
child2(added)   = kc;
child2(removed) = k0;
Cnow = interior & (child2==kc);
n_after_edits = nnz(Cnow);

%% small residual correction to hit N exactly (bounded: only the count mismatch, not a redesign)
if n_after_edits ~= N
    Nm = interior & (child2==k0);
    phi = signed_distance_field(Cnow, Nm, max_dist, 'octagon');
    D2 = interior & (child2==k0 | child2==kc);
    idx = find(D2);
    [~,ord] = sort(phi(idx) + 1e-3*rand(numel(idx),1), 'descend');
    sel = false(numel(idx),1); sel(ord(1:min(N,numel(idx)))) = true;
    child2(idx) = k0;
    child2(idx(sel)) = kc;
end

if nargout > 1
    info = struct('count_in',nC,'count_after_edits',n_after_edits,'count_out',nnz(child2==kc & interior), ...
                  'residual', n_after_edits - N, 'cells_added', nnz(added), 'cells_removed', nnz(removed));
end
end

%% ============================================================================================
%% Local functions below -- private to this file, nothing else on the path is required.
%% ============================================================================================

function S = skeletonize(BW)
% Zhang-Suen thinning, vectorised, base MATLAB/Octave, no Image Processing Toolbox needed.
% BW: logical matrix. Returns S: logical, 1-pixel-wide skeleton, same size.
[h,w] = size(BW);
I = false(h+2,w+2);
I(2:h+1,2:w+1) = logical(BW);
changed = true;
iter = 0;
while changed && iter < 200
    changed = false; iter = iter+1;
    for sub = 1:2
        P2 = I([1 1:end-1],:);              % N
        P6 = I([2:end end],:);              % S
        P4 = I(:,[2:end end]);              % E
        P8 = I(:,[1 1:end-1]);              % W
        P3 = I([1 1:end-1],[2:end end]);    % NE
        P5 = I([2:end end],[2:end end]);    % SE
        P7 = I([2:end end],[1 1:end-1]);    % SW
        P9 = I([1 1:end-1],[1 1:end-1]);    % NW
        B = double(P2)+P3+P4+P5+P6+P7+P8+P9;
        seq = {P2,P3,P4,P5,P6,P7,P8,P9,P2};
        A = zeros(size(I));
        for k = 1:8
            A = A + (~seq{k} & seq{k+1});
        end
        if sub==1
            cond = (~(P2 & P4 & P6)) & (~(P4 & P6 & P8));
        else
            cond = (~(P2 & P4 & P8)) & (~(P2 & P6 & P8));
        end
        marked = I & (B>=2) & (B<=6) & (A==1) & cond;
        if any(marked(:))
            I(marked) = false;
            changed = true;
        end
    end
end
S = I(2:h+1,2:w+1);

%% safety net: parallel (simultaneous) deletion can occasionally erase an entire small component
%% in one pass (e.g. every pixel of a 2x2 block can satisfy the deletion test at once). Restore one
%% representative pixel -- the most interior one -- for any component that lost all its pixels, so
%% no material silently disappears.
L = label_components(BW, 8);
if max(L(:)) > 0
    Dsafe = local_thickness(BW);        % distance to background, used to pick the "most interior" pixel
    for lbl = 1:max(L(:))
        comp = (L==lbl);
        if ~any(S(comp))
            Dc = Dsafe; Dc(~comp) = -Inf;
            [~,imax] = max(Dc(:));
            S(imax) = true;
        end
    end
end
end

function L = label_components(mask, conn)
% Connected-component labeling, no toolbox needed, vectorised (min-propagation: every foreground
% pixel repeatedly takes the minimum label among itself and its foreground neighbours, until no
% pixel changes; each component converges to the smallest linear index inside it). conn = 4 or 8
% (default 8). Labels returned are 0 (background) or 1..k in no particular order.
if nargin < 2; conn = 8; end
[h,w] = size(mask);
mask = logical(mask);
idx = reshape(1:h*w, h, w);
L = zeros(h,w); L(mask) = idx(mask);
if conn==4
    offs = [-1 0; 1 0; 0 -1; 0 1];
else
    offs = [-1 0; 1 0; 0 -1; 0 1; -1 -1; -1 1; 1 -1; 1 1];
end
changed = true; iter = 0; max_iter = 4*(h+w);
while changed && iter < max_iter
    changed = false; iter = iter+1;
    for k = 1:size(offs,1)
        dr = offs(k,1); dc = offs(k,2);
        shifted = zeros(h,w);
        rs = max(1,1-dr):min(h,h-dr); cs = max(1,1-dc):min(w,w-dc);
        shifted(rs,cs) = L(rs+dr,cs+dc);
        cand = mask & shifted>0 & (shifted<L | L==0);
        if any(cand(:))
            L(cand) = shifted(cand);
            changed = true;
        end
    end
end
[~,~,ic] = unique(L(mask));
Lc = zeros(h,w); Lc(mask) = ic;
L = Lc;
end

function D = local_thickness(BW)
% distance (chessboard-ish, 4/8 alternating) from every foreground pixel to the nearest background
% pixel, 0 outside BW. No toolbox needed.
[h,w] = size(BW);
K4 = [0 1 0; 1 1 1; 0 1 0]; K8 = ones(3);
bg = ~BW;
reach = bg; D = zeros(h,w); maxd = max(h,w);
for k = 1:maxd
    if mod(k,2)==1; K = K4; else; K = K8; end
    grown = conv2(double(reach), K, 'same') > 0;
    front = grown & ~reach;
    if ~any(front(:)); break; end
    D(front) = k;
    reach = grown;
end
D(bg) = 0;
end

function [ep, bp, deg] = skeleton_endpoints(S)
% ep: linear indices of endpoint pixels (degree 1, 8-connectivity)
% bp: linear indices of branch-point pixels (degree >= 3)
% deg: degree map (same size as S)
[h,w] = size(S);
nb = zeros(h,w);
for dr = -1:1
    for dc = -1:1
        if dr==0 && dc==0; continue; end
        shifted = false(h,w);
        rs = max(1,1-dr):min(h,h-dr); cs = max(1,1-dc):min(w,w-dc);
        shifted(rs,cs) = S(rs+dr,cs+dc);
        nb = nb + double(shifted);
    end
end
deg = nb; deg(~S) = 0;
ep = find(S & deg==1);
bp = find(S & deg>=3);
end

function path = walk_skeleton(S, start_idx, max_steps)
% Walk a simple path on the skeleton S starting at start_idx (must be an endpoint or any degree<=2
% pixel), choosing at each step the unique unvisited 8-neighbour still on the skeleton. Stops when
% no unvisited neighbour remains (reached the other end / a branch point) or after max_steps.
% path: column vector of linear indices, path(1) = start_idx.
[h,w] = size(S);
[r0,c0] = ind2sub([h w], start_idx);
path = start_idx;
visited = false(h,w); visited(start_idx) = true;
r = r0; c = c0;
for step = 1:max_steps
    found = false;
    for dr = -1:1
        for dc = -1:1
            if dr==0 && dc==0; continue; end
            r1 = r+dr; c1 = c+dc;
            if r1>=1 && r1<=h && c1>=1 && c1<=w && S(r1,c1) && ~visited(r1,c1)
                r = r1; c = c1; visited(r,c) = true;
                path(end+1,1) = sub2ind([h w], r, c); %#ok<AGROW>
                found = true; break;
            end
        end
        if found; break; end
    end
    if ~found; break; end
    % stop AFTER including a branch point (degree>=3 among 8-neighbours on S), so the caller can see it
    degk = 0;
    for dr=-1:1, for dc=-1:1
        if dr==0 && dc==0; continue; end
        r1=r+dr; c1=c+dc;
        if r1>=1 && r1<=h && c1>=1 && c1<=w && S(r1,c1); degk = degk+1; end
    end, end
    if degk >= 3; break; end
end
end

function phi = signed_distance_field(Cmask, Nmask, max_dist, metric)
% phi > 0 inside conductive matter, phi < 0 in non-conductive cells; |phi| = 0.5 on the two cell
% layers touching the interface, 1.5 on the next ones, etc. Cells in neither mask get phi = 0.
if nargin < 3 || isempty(max_dist); max_dist = 12; end
if nargin < 4 || isempty(metric);   metric = 'octagon'; end
d_to_C = dist_to(Cmask, max_dist, metric);   % outside cells: distance to nearest conductive cell
d_to_N = dist_to(Nmask, max_dist, metric);   % inside cells : distance to nearest non-conductive cell
phi = zeros(size(Cmask));
phi(Cmask) =   d_to_N(Cmask) - 0.5;
phi(Nmask) = -(d_to_C(Nmask) - 0.5);
end

function d = dist_to(src, max_dist, metric)
% distance (in cells) from every cell to the nearest true cell of src, 0 on src, saturated at
% max_dist+1. metric = 'euclid' uses bwdist (Image Processing Toolbox); default 'octagon' needs
% no toolbox.
if strcmp(metric,'euclid')
    d = min(bwdist(src), max_dist+1);
    return
end
K4 = [0 1 0; 1 1 1; 0 1 0];
K8 = ones(3);
reach = src;
d = (max_dist+1)*ones(size(src));
d(src) = 0;
for k = 1:max_dist
    if mod(k,2) == 1; K = K4; else K = K8; end
    grown = conv2(double(reach), K, 'same') > 0;
    front = grown & ~reach;
    if ~any(front(:)); break; end
    d(front) = k;
    reach = grown;
end
end

function [dr,dc] = local_direction(path, h, w)
% outward direction estimate from a short walked path (path(1) = the endpoint)
if numel(path) < 2; dr = NaN; dc = NaN; return; end
[r1,c1] = ind2sub([h w], path(1));
k = min(numel(path), 5);
[r2,c2] = ind2sub([h w], path(k));
dr = r1-r2; dc = c1-c2;
if dr==0 && dc==0; dr = NaN; dc = NaN; end
end

function [dr,dc] = local_tangent(S, idx, h, w)
% local tangent direction at a skeleton pixel: average vector to its skeleton neighbours
[r0,c0] = ind2sub([h w], idx);
dr = 0; dc = 0; n = 0;
for a = -1:1
    for b = -1:1
        if a==0 && b==0; continue; end
        r1=r0+a; c1=c0+b;
        if r1>=1 && r1<=h && c1>=1 && c1<=w && S(r1,c1); dr=dr+a; dc=dc+b; n=n+1; end
    end
end
if n==0; dr = randn; dc = randn; end
end

function [mask,changed,S] = draw_segment(mask, changed, S, start_idx, ang, L, radius, h, w, interior)
[r0,c0] = ind2sub([h w], start_idx);
disk = octagon_disk_local(max(1,round(radius)));
for t = 1:L
    r = round(r0 + t*sin(ang)); c = round(c0 + t*cos(ang));
    if r<1 || r>h || c<1 || c>w; break; end
    p = sub2ind([h w], r, c);
    [mask,changed] = stamp(mask, changed, p, disk, h, w, interior, true);
    S(p) = true;                                     % the drawn centreline is itself already thin
end
end

function [mask,changed] = stamp(mask, changed, center_idx, disk, h, w, interior, value)
% write 'value' into mask under 'disk' centred at center_idx, restricted to interior design cells
[r0,c0] = ind2sub([h w], center_idx);
r = size(disk,1); rad = (r-1)/2;
rs = max(1,r0-rad):min(h,r0+rad); cs = max(1,c0-rad):min(w,c0+rad);
drs = rs - (r0-rad) + 1; dcs = cs - (c0-rad) + 1;
sub = disk(drs,dcs) & interior(rs,cs);
region = mask(rs,cs);
if value
    changed(rs,cs) = changed(rs,cs) | (sub & ~region);
    region(sub) = true;
else
    changed(rs,cs) = changed(rs,cs) | (sub & region);
    region(sub) = false;
end
mask(rs,cs) = region;
end

function disk = octagon_disk_local(r)
persistent cache
if isempty(cache); cache = containers.Map('KeyType','double','ValueType','any'); end
if isKey(cache, r); disk = cache(r); return; end
n = 2*r+1; P = false(n,n); P(r+1,r+1) = true;
K4 = [0 1 0; 1 1 1; 0 1 0]; K8 = ones(3); reach = P;
for k = 1:r
    if mod(k,2)==1; K = K4; else; K = K8; end
    reach = conv2(double(reach), K, 'same') > 0;
end
disk = reach;
cache(r) = disk;
end

function v = get_opt(opts, name, default)
if isfield(opts, name) && ~isempty(opts.(name)); v = opts.(name); else; v = default; end
end
