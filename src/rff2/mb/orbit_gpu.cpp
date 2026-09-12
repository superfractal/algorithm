// Created by GPT-6 on 2026-09-11
// Modified by GPT-6 on 2026-09-12
#include "orbit_gpu.hpp"
#include <bit>
#include <filesystem>
#include <cstdlib>
unsigned OrbitGPU::configuredFlags() {
    const char *v = std::getenv("RFF_GPU_EXPERIMENT");
    return v ? std::stoul(v) : 15u;
}
std::string OrbitGPU::shaderFile(const char *p, unsigned f) {
    return (std::filesystem::path(p).parent_path() / ("orbit_" + std::to_string(f & 7) + ".spv"))
        .string();
}
void OrbitGPU::orbitPlan(unsigned bits, unsigned steps, unsigned count, bool first, bool last) {
    if (fraction == bits && iterations == steps && segments == count && uploadState == first &&
        downloadState == last)
        return;
    if (bits < 32 || bits > 16000000 || !steps || steps > 32 || !count)
        throw std::runtime_error("resident orbit range");
    unsigned n = std::max(512u, std::bit_ceil(2 * ((bits + 47) / 16))), jobs = products() * count;
    if (4ull * (n / 2) * 65535 * 65535 >= uint64_t(moduli[0]) * moduli[1])
        throw std::runtime_error("signed CRT coefficient bound");
    VkDeviceSize meta = (products() + 3ull) * count * n * 4;
    if (fraction != bits || segments != count) {
        if (!first)
            throw std::runtime_error("resident shape change");
        resetPlan();
        plan(n, jobs);
        fraction = iterations = segments = 0;
        free(bufs[4]);
        free(download);
        bufs[4] = {};
        download = {};
        VkDeviceSize bytes = meta + (12ull * count + 20ull * 32 * count) * 4;
        if (bytes > props.limits.maxStorageBufferRange)
            throw std::runtime_error("resident storage range");
        bufs[4] = buffer(bytes, false);
        download = buffer(bytes, true);
        if (upload.size < bytes)
            throw std::runtime_error("resident upload capacity");
        VkDescriptorBufferInfo info{bufs[4].b, 0, bytes};
        VkWriteDescriptorSet w{VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET};
        w.dstSet = descriptor;
        w.dstBinding = 4;
        w.descriptorCount = 1;
        w.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        w.pBufferInfo = &info;
        vkUpdateDescriptorSets(device, 1, &w, 0, nullptr);
    }
    phases.clear();
    ck(vkResetCommandBuffer(command, 0));
    VkCommandBufferBeginInfo begin{VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    ck(vkBeginCommandBuffer(command, &begin));
    vkCmdResetQueryPool(command, queries, 0, profile ? 4096 : 2);
    if (first) {
        VkBufferCopy in[2] = {
            {jobs * uint64_t(n) * 4, jobs * uint64_t(n) * 4, 2ull * count * n * 4},
            {meta, meta, 12ull * count * 4}};
        vkCmdCopyBuffer(command, upload.b, bufs[4].b, 2, in);
    }
    barrier(command,
            VK_ACCESS_TRANSFER_WRITE_BIT | VK_ACCESS_SHADER_WRITE_BIT,
            VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT,
            VK_PIPELINE_STAGE_TRANSFER_BIT | VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
            VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT);
    vkCmdBindPipeline(command, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline);
    vkCmdBindDescriptorSets(
        command, VK_PIPELINE_BIND_POINT_COMPUTE, layout, 0, 1, &descriptor, 0, nullptr);
    vkCmdWriteTimestamp(command, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, queries, 0);
    uint32_t inv[2];
    for (int p = 0; p < 2; ++p)
        inv[p] =
            uint64_t(power(n, moduli[p] - 2, moduli[p])) * ((1ull << 32) % moduli[p]) % moduli[p];
    auto dispatch = [&](unsigned phase,
                        uint64_t tasks,
                        unsigned span = 0,
                        unsigned offset = 0,
                        unsigned entries = 0,
                        unsigned overrideJobs = 0) {
        uint32_t push[] = {
            n, overrideJobs ? overrideJobs : jobs, phase, span, inv[0], inv[1], offset, entries};
        vkCmdPushConstants(command, layout, VK_SHADER_STAGE_COMPUTE_BIT, 0, 32, push);
        unsigned groups = (tasks + 255) / 256;
        vkCmdDispatch(command, std::min(65535u, groups), (groups + 65534) / 65535, 1);
        barrier(command,
                VK_ACCESS_SHADER_WRITE_BIT,
                VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT,
                VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT);
        if (profile) {
            phases.push_back(phase);
            if (phases.size() + 2 >= 4096)
                throw std::runtime_error("profile capacity");
            vkCmdWriteTimestamp(
                command, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, queries, 1 + phases.size());
        }
    };
    auto carry = [&](unsigned len, unsigned cnt, bool sign = false) {
        unsigned blocks = len / 256;
        dispatch(sign ? 32 : 8, cnt * uint64_t((blocks + 255) / 256) * 256, 0, len, blocks, cnt);
        if (blocks > 256)
            dispatch(sign ? 32 : 8, cnt * 256ull, 0, len + blocks, (blocks + 255) / 256, cnt);
    };
    for (unsigned step = 0; step < steps; ++step) {
        dispatch(19, count);
        dispatch(20, count * uint64_t(n));
        dispatch(28, count * 2ull, step);
        for (unsigned len = n; len > 512;) {
            bool fuse = (flags & 4) && len >= 2048;
            dispatch(fuse ? 40 : 1, count * uint64_t(n) * (fuse ? 1 : 2), len);
            len /= fuse ? 4 : 2;
        }
        dispatch(5, count * 2ull * n);
        dispatch(2, count * 2ull * n);
        dispatch(6, jobs * uint64_t(n));
        for (unsigned len = 1024; len <= n;) {
            bool fuse = (flags & 4) && len <= n / 2;
            dispatch(fuse ? 41 : 3, jobs * uint64_t(n) / (fuse ? 2 : 1), fuse ? len * 2 : len);
            len *= fuse ? 4 : 2;
        }
        dispatch(4, jobs * uint64_t(n));
        dispatch((flags & 2) ? 31 : 7, jobs * uint64_t(n));
        carry(n, jobs, (flags & 2) != 0);
        dispatch((flags & 2) ? 33 : 9, jobs * uint64_t(n));
        if (!(flags & 2)) {
            dispatch(21, count * uint64_t(n));
            carry(n, count);
            dispatch(22, count * uint64_t(n));
        }
        dispatch(23, count * uint64_t(n), 0, bits);
        dispatch(26, count * 2ull);
        dispatch(24, count * uint64_t(n));
        carry(n / 2, count * 2);
        dispatch(25, count * uint64_t(n));
        dispatch(27, count * uint64_t(n));
    }
    vkCmdWriteTimestamp(command, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, queries, 1);
    barrier(command,
            VK_ACCESS_SHADER_WRITE_BIT,
            VK_ACCESS_TRANSFER_READ_BIT,
            VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
            VK_PIPELINE_STAGE_TRANSFER_BIT);
    VkBufferCopy out[2] = {
        {jobs * uint64_t(n) * 4, jobs * uint64_t(n) * 4, count * uint64_t(n) * 4},
        {meta, meta, (12ull * count + 20ull * steps * count) * 4}};
    vkCmdCopyBuffer(command, bufs[4].b, download.b, last ? 2 : 1, last ? out : out + 1);
    barrier(command,
            VK_ACCESS_TRANSFER_WRITE_BIT,
            VK_ACCESS_HOST_READ_BIT,
            VK_PIPELINE_STAGE_TRANSFER_BIT,
            VK_PIPELINE_STAGE_HOST_BIT);
    ck(vkEndCommandBuffer(command));
    fraction = bits;
    iterations = steps;
    segments = count;
    uploadState = first;
    downloadState = last;
}
void OrbitGPU::submitOrbit(unsigned bits,
                           unsigned steps,
                           const std::vector<Segment> &jobs,
                           bool first,
                           bool last,
                           const std::function<bool()> &cancel) {
    if (pending)
        throw std::runtime_error("orbit submission already pending");
    if (!first && !resident)
        throw std::runtime_error("no resident session");
    callStart = std::chrono::steady_clock::now();
    if (cancel()) {
        resident = false;
        throw std::runtime_error("GPU orbit cancelled");
    }
    orbitPlan(bits, steps, jobs.size(), first, last);
    size_t n = planSize, L = n / 2, S = jobs.size(), st = products() * S * n, cs = st + 2 * S * L,
           meta = cs + 4 * S * L;
    if (first) {
        auto *in = static_cast<uint32_t *>(upload.mapped);
        std::memset(in, 0, bufs[4].size);
        std::vector<uint16_t> digits(L);
        for (size_t j = 0; j < S; ++j) {
            mpz_srcptr v[] = {jobs[j].x, jobs[j].y, jobs[j].cr, jobs[j].ci};
            for (unsigned c = 0; c < 4; ++c) {
                if (mpz_sizeinbase(v[c], 2) > bits + 16)
                    throw std::runtime_error("orbit input integer range");
                size_t written = 0;
                mpz_export(digits.data(), &written, -1, 2, 0, 0, v[c]);
                auto *dst = in + (c < 2 ? st : cs) + (j * 2 + c % 2) * L;
                for (size_t i = 0; i < written; ++i)
                    dst[i] = digits[i];
                in[meta + 12 * j + c] = mpz_sgn(v[c]) < 0;
            }
        }
    }
    if (cancel()) {
        resident = false;
        throw std::runtime_error("GPU orbit cancelled");
    }
    pendingJobs = jobs;
    VkSubmitInfo submit{VK_STRUCTURE_TYPE_SUBMIT_INFO};
    submit.commandBufferCount = 1;
    submit.pCommandBuffers = &command;
    ck(vkResetFences(device, 1, &fence));
    submitStart = std::chrono::steady_clock::now();
    ck(vkQueueSubmit(queue, 1, &submit, fence));
    pending = true;
    resident = !last;
}
NTTGPU::Timing OrbitGPU::finishOrbit(const std::function<bool()> &cancel) {
    if (!pending)
        throw std::runtime_error("no pending orbit");
    size_t n = planSize, L = n / 2, S = pendingJobs.size(), st = products() * S * n,
           meta = st + 3 * S * n;
    bool stopped = false;
    for (;;) {
        auto r = vkWaitForFences(device, 1, &fence, 1, 10000000);
        stopped = stopped || cancel();
        if (r == VK_SUCCESS)
            break;
        if (r != VK_TIMEOUT)
            ck(r);
    }
    pending = false;
    Timing t;
    t.submit =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - submitStart).count();
    if (stopped) {
        resident = false;
        throw std::runtime_error("GPU orbit cancelled");
    }
    uint64_t ticks[2];
    ck(vkGetQueryPoolResults(device,
                             queries,
                             0,
                             2,
                             sizeof(ticks),
                             ticks,
                             8,
                             VK_QUERY_RESULT_64_BIT | VK_QUERY_RESULT_WAIT_BIT));
    t.gpu = (ticks[1] - ticks[0]) * props.limits.timestampPeriod * 1e-9;
    phaseSeconds.assign(64, 0);
    if (profile) {
        std::vector<uint64_t> marks(phases.size());
        ck(vkGetQueryPoolResults(device,
                                 queries,
                                 2,
                                 marks.size(),
                                 marks.size() * 8,
                                 marks.data(),
                                 8,
                                 VK_QUERY_RESULT_64_BIT | VK_QUERY_RESULT_WAIT_BIT));
        uint64_t prev = ticks[0];
        for (size_t k = 0; k < marks.size(); ++k) {
            phaseSeconds[phases[k]] += (marks[k] - prev) * props.limits.timestampPeriod * 1e-9;
            prev = marks[k];
        }
    }
    auto *out = static_cast<uint32_t *>(download.mapped);
    for (size_t j = 0; j < S; ++j)
        if (out[meta + 12 * j + 7]) {
            resident = false;
            throw std::runtime_error("GPU orbit integer overflow");
        }
    trace.assign(out + meta + 12 * S, out + meta + 12 * S + 20ull * iterations * S);
    if (downloadState) {
        std::vector<uint16_t> digits(L);
        for (size_t j = 0; j < S; ++j)
            for (unsigned c = 0; c < 2; ++c) {
                for (size_t i = 0; i < L; ++i)
                    digits[i] = uint16_t(out[st + (2 * j + c) * L + i]);
                auto target = c ? pendingJobs[j].resultY : pendingJobs[j].resultX;
                mpz_import(target, L, -1, 2, 0, 0, digits.data());
                if (out[meta + 12 * j + c])
                    mpz_neg(target, target);
            }
    }
    t.total = std::chrono::duration<double>(std::chrono::steady_clock::now() - callStart).count();
    return t;
}
NTTGPU::Timing OrbitGPU::orbit(unsigned b,
                               unsigned s,
                               const std::vector<Segment> &j,
                               const std::function<bool()> &c) {
    submitOrbit(b, s, j, true, true, c);
    return finishOrbit(c);
}
