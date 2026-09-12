// Created by GPT-6 on 2026-09-12
// Modified by GPT-6 on 2026-09-12
// Phases 40/41 fuse two transform stages; preserve butterfly and carry composition order.
if (p.phase == 40 || p.phase == 41) {
    bool inv = p.phase == 41;
    uint streams = inv ? p.jobs * 2 : (p.jobs / PRODUCTS) * 4;
    if (t >= streams * (p.n / 4))
        return;
    uint stream = t / (p.n / 4), ix = t % (p.n / 4), prime = inv ? stream % 2 : (stream % 4) / 2;
    uint base = (inv ? stream * 2
                     : ((stream / 4) * PRODUCTS * 4 + (stream % 4) / 2 * 2 + (stream % 2) * 4)) *
                p.n;
    uint q = p.span / 4, j = ix % q, at = base + (ix / q) * p.span + j,
         ro = prime * p.n + (inv ? p.n / 2 : 0);
    uint w0 = roots[ro + j * (p.n / p.span)], w1 = roots[ro + (j + q) * (p.n / p.span)],
         w2 = roots[ro + j * (p.n / (p.span / 2))];
    uint a = data[at], b = data[at + q], c = data[at + 2 * q], d = data[at + 3 * q];
    uint u, v, x, y;
    if (!inv) {
        u = plus(a, c, prime);
        v = plus(b, d, prime);
        x = mont(minus(a, c, prime), w0, prime);
        y = mont(minus(b, d, prime), w1, prime);
        data[at] = plus(u, v, prime);
        data[at + q] = mont(minus(u, v, prime), w2, prime);
        data[at + 2 * q] = plus(x, y, prime);
        data[at + 3 * q] = mont(minus(x, y, prime), w2, prime);
    } else {
        b = mont(b, w2, prime);
        d = mont(d, w2, prime);
        u = plus(a, b, prime);
        v = minus(a, b, prime);
        x = mont(plus(c, d, prime), w0, prime);
        y = mont(minus(c, d, prime), w1, prime);
        data[at] = plus(u, x, prime);
        data[at + q] = plus(v, y, prime);
        data[at + 2 * q] = minus(u, x, prime);
        data[at + 3 * q] = minus(v, y, prime);
    }
    return;
}
if (p.phase == 31) {
    if (t >= p.jobs * p.n)
        return;
    uint job = t / p.n, i = t % p.n;
    int64_t sum = int64_t(coefficients[t] & 65535ul);
    if (i >= 1)
        sum += int64_t((coefficients[t - 1] >> 16) & 65535ul);
    if (i >= 2)
        sum += int64_t((coefficients[t - 2] >> 32) & 65535ul);
    if (i >= 3)
        sum += int64_t(coefficients[t - 3]) >> 48;
    result[t] = uint(sum);
    uint f = 0;
    for (uint k = 0; k < 5; ++k)
        f |= uint(((sum + int64_t(k) - 1) >> 16) + 1) << (3 * k);
    tile[lane] = f;
    barrier();
    signedScan(lane);
    uint base = job * mapStride();
    maps[base + i] = tile[lane];
    if (lane == 255)
        maps[base + p.n + i / 256] = tile[lane];
    return;
}
if (p.phase == 32) {
    uint groups = (p.count + 255) / 256;
    if (group >= p.jobs * groups)
        return;
    uint job = group / groups, g = group % groups, i = g * 256 + lane,
         base = job * mapStride() + p.offset;
    tile[lane] = i < p.count ? maps[base + i] : 18056u;
    barrier();
    signedScan(lane);
    if (i < p.count)
        maps[base + i] = tile[lane];
    if (lane == 255)
        maps[base + p.count + g] = tile[lane];
    return;
}
if (p.phase == 33) {
    if (t >= p.jobs * p.n)
        return;
    uint job = t / p.n, i = t % p.n, block = i / 256, grand = block / 256, base = job * mapStride(),
         blocks = p.n / 256, c = 1;
    if (grand > 0)
        c = signedCarry(maps[base + p.n + blocks + grand - 1], c);
    if (block % 256 > 0)
        c = signedCarry(maps[base + p.n + block - 1], c);
    if (i % 256 > 0)
        c = signedCarry(maps[base + i - 1], c);
    result[t] = (result[t] + uint(int(c) - 1)) & 65535u;
    return;
}
