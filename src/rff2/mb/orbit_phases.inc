// Created by GPT-6 on 2026-09-11
// Modified by GPT-6 on 2026-09-12
// Result layout: products[J*N], state[2*S*L], c[2*S*L], temporary[2*S*L], metadata[12*S].
// Each numbered phase is a separate host dispatch with a memory barrier between phases.
if (p.phase >= 19) {
    uint segments = p.jobs / 3, L = p.n / 2, st = p.jobs * p.n, cs = st + 2 * segments * L,
         tmp = cs + 2 * segments * L, meta = tmp + 2 * segments * L;
    if (p.phase == 19) {
        if (t >= segments)
            return;
        result[meta + 12 * t + 4] = 0;
        result[meta + 12 * t + 5] = 0;
        result[meta + 12 * t + 6] = 0;
        result[meta + 12 * t + 8] = 0;
        result[meta + 12 * t + 9] = 0;
        return;
    }
    if (p.phase == 20) {
        if (t >= segments * p.n)
            return;
        uint seg = t / p.n, i = t % p.n;
        uint x = i < L ? result[st + 2 * seg * L + i] : 0,
             y = i < L ? result[st + (2 * seg + 1) * L + i] : 0;
        if (x != y)
            atomicMax(result[meta + 12 * seg + 4], i + 1);
        if (x != 0)
            atomicMax(result[meta + 12 * seg + 8], i + 1);
        if (y != 0)
            atomicMax(result[meta + 12 * seg + 9], i + 1);
        for (uint prime = 0; prime < 2; ++prime) {
            uint r2 = prime == 0 ? 1172168163u : 959408210u;
            data[(seg * 12 + prime * 2) * p.n + i] = mont(x, r2, prime);
            data[(seg * 12 + 4 + prime * 2) * p.n + i] = mont(y, r2, prime);
        }
        return;
    }
    if (p.phase == 21) {
        if (t >= segments * p.n)
            return;
        uint seg = t / p.n, i = t % p.n, at = result[meta + 12 * seg + 4];
        bool realNegative =
            at != 0 && result[st + 2 * seg * L + at - 1] < result[st + (2 * seg + 1) * L + at - 1];
        uint a = result[(3 * seg + uint(realNegative)) * p.n + i],
             b = result[(3 * seg + uint(!realNegative)) * p.n + i];
        result[3 * seg * p.n + i] = uint(int(a) - int(b));
        uint f = uint(a < b) | (uint(a <= b) << 2) | (uint(a <= b) << 4) | (uint(a <= b) << 6);
        tile[lane] = f;
        barrier();
        scan(lane);
        uint base = seg * mapStride();
        maps[base + i] = tile[lane];
        if (lane == 255)
            maps[base + p.n + i / 256] = tile[lane];
        return;
    }
    if (p.phase == 22) {
        if (t >= segments * p.n)
            return;
        uint seg = t / p.n, i = t % p.n, block = i / 256, grand = block / 256,
             base = seg * mapStride(), blocks = p.n / 256, c = 0;
        if (grand > 0)
            c = applyCarry(maps[base + p.n + blocks + grand - 1], c);
        if (block % 256 > 0)
            c = applyCarry(maps[base + p.n + block - 1], c);
        if (i % 256 > 0)
            c = applyCarry(maps[base + i - 1], c);
        result[3 * seg * p.n + i] = (result[3 * seg * p.n + i] - c) & 65535u;
        return;
    }
    // Shift the exact integer products to the fixed-point orbit scale.
    if (p.phase == 23) {
        if (t >= segments * 2 * L)
            return;
        uint comp = t / L, seg = comp / 2, c = comp % 2, i = t % L;
        uint shift = p.offset - c, source = (3 * seg + (c == 0 ? 0 : 2)) * p.n,
             word = i + shift / 16, r = shift % 16;
        uint v = word < p.n ? result[source + word] >> r : 0;
        if (r != 0 && word + 1 < p.n)
            v |= (result[source + word + 1] << (16 - r)) & 65535u;
        result[tmp + t] = v;
        if (v != result[cs + t])
            atomicMax(result[meta + 12 * seg + 5 + c], i + 1);
        word += L;
        uint lost = word < p.n ? result[source + word] >> r : 0;
        if (r != 0 && word + 1 < p.n)
            lost |= (result[source + word + 1] << (16 - r)) & 65535u;
        if (lost != 0)
            atomicOr(result[meta + 12 * seg + 7], 1u);
        return;
    }
    if (p.phase == 24) {
        if (t >= segments * 2 * L)
            return;
        uint comp = t / L, seg = comp / 2, c = comp % 2, i = t % L;
        uint flags = result[meta + 12 * seg + 10 + c];
        bool subtract = (flags & 1) != 0, swapOrder = (flags & 2) != 0;
        uint a = result[(swapOrder ? cs : tmp) + t], b = result[(swapOrder ? tmp : cs) + t];
        uint v = subtract ? uint(int(a) - int(b)) : a + b, f = 0;
        for (uint k = 0; k < 4; ++k) {
            uint carry = subtract ? uint(a < b + min(k, 1u)) : ((v + min(k, 1u)) >> 16);
            f |= carry << (2 * k);
        }
        result[tmp + t] = v;
        tile[lane] = f;
        barrier();
        scan(lane);
        uint base = comp * mapStride();
        maps[base + i] = tile[lane];
        if (lane == 255)
            maps[base + L + i / 256] = tile[lane];
        // Signs must remain immutable while all workgroups read the old orbit.
        return;
    }
    if (p.phase == 25) {
        if (t >= segments * 2 * L)
            return;
        uint comp = t / L, seg = comp / 2, c = comp % 2, i = t % L;
        bool subtract = (result[meta + 12 * seg + 10 + c] & 1) != 0;
        uint block = i / 256, grand = block / 256, base = comp * mapStride(), blocks = L / 256,
             carry = 0;
        if (grand > 0)
            carry = applyCarry(maps[base + L + blocks + grand - 1], carry);
        if (block % 256 > 0)
            carry = applyCarry(maps[base + L + block - 1], carry);
        if (i % 256 > 0)
            carry = applyCarry(maps[base + i - 1], carry);
        result[tmp + t] = (result[tmp + t] + (subtract ? uint(-int(carry)) : carry)) & 65535u;
        return;
    }
    if (p.phase == 26) {
        if (t >= segments * 2)
            return;
        uint comp = t, seg = comp / 2, c = comp % 2, xyAt = result[meta + 12 * seg + 4];
        bool realNegative = xyAt != 0 && result[st + 2 * seg * L + xyAt - 1] <
                                             result[st + (2 * seg + 1) * L + xyAt - 1];
        uint sg =
            c == 0 ? uint(realNegative) : (result[meta + 12 * seg] ^ result[meta + 12 * seg + 1]);
        bool subtract = sg != result[meta + 12 * seg + 2 + c];
        uint cmp = result[meta + 12 * seg + 5 + c];
        bool swapOrder = subtract && cmp != 0 &&
                         result[tmp + comp * L + cmp - 1] < result[cs + comp * L + cmp - 1];
        result[meta + 12 * seg + 8 + c] = swapOrder ? result[meta + 12 * seg + 2 + c] : sg;
        result[meta + 12 * seg + 10 + c] = uint(subtract) | uint(swapOrder) * 2;
        return;
    }
    if (p.phase == 27) {
        if (t >= segments * 2 * L)
            return;
        uint comp = t / L, seg = comp / 2, c = comp % 2, i = t % L;
        result[st + t] = result[tmp + t];
        if (i == 0)
            result[meta + 12 * seg + c] = result[meta + 12 * seg + 8 + c];
        return;
    }
    // Export the leading eight limbs and discarded-limb count for CPU error bounds.
    if (p.phase == 28) {
        if (t >= segments * 2)
            return;
        uint comp = t, seg = comp / 2, c = comp % 2, msd = result[meta + 12 * seg + 8 + c];
        uint cut = msd > 8 ? msd - 8 : 0,
             at = meta + 12 * segments + (p.span * segments * 2 + comp) * 10;
        for (uint k = 0; k < 8; ++k)
            result[at + k] = cut + k < L ? result[st + comp * L + cut + k] : 0;
        result[at + 8] = cut;
        result[at + 9] = msd != 0 ? result[meta + 12 * seg + c] : 0u;
        return;
    }
}
