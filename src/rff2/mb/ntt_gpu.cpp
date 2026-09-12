// Created by GPT-6 on 2026-09-11
// Modified by GPT-6 on 2026-09-12
#include "ntt_gpu.hpp"
#include <bit>
void NTTGPU::resetPlan() {
    if (fence)
        vkDestroyFence(device, fence, nullptr);
    fence = {};
    if (queries)
        vkDestroyQueryPool(device, queries, nullptr);
    queries = {};
    if (command)
        vkFreeCommandBuffers(device, pool, 1, &command);
    command = {};
    for (auto &b: bufs)
        if (b.b) {
            free(b);
            b = {};
        }
    for (auto *b: {&upload, &download, &rootUpload})
        if (b->b) {
            free(*b);
            *b = {};
        }
    planSize = planJobs = 0;
}
void NTTGPU::plan(unsigned n, unsigned jobs) {
    if (n == planSize && jobs == planJobs)
        return;
    if (!jobs || uint64_t(jobs) * n * 52 + uint64_t(n) * 16 > 1024ull * 1024 * 1024)
        throw std::runtime_error("NTT memory budget requires smaller job tile");
    resetPlan();
    prepare(n);
    const size_t blocks = n / 256, stride = n + blocks + (blocks + 255) / 256 + 1;
    const VkDeviceSize bytes[] = {jobs * 4ull * n * 4,
                                  2ull * n * 4,
                                  jobs * uint64_t(n) * 8,
                                  jobs * stride * 4,
                                  jobs * uint64_t(n) * 4};
    for (unsigned i = 0; i < 5; ++i) {
        if (bytes[i] > props.limits.maxStorageBufferRange)
            throw std::runtime_error("storage range requires smaller tile");
        bufs[i] = buffer(bytes[i], false);
    }
    upload = buffer(bytes[0], true);
    download = buffer(bytes[4], true);
    rootUpload = buffer(bytes[1], true);
    std::memcpy(rootUpload.mapped, rootTable.data(), bytes[1]);
    ck(vkResetDescriptorPool(device, descriptors, 0));
    VkDescriptorSetAllocateInfo da{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
    da.descriptorPool = descriptors;
    da.descriptorSetCount = 1;
    da.pSetLayouts = &setLayout;
    ck(vkAllocateDescriptorSets(device, &da, &descriptor));
    VkDescriptorBufferInfo infos[5];
    VkWriteDescriptorSet writes[5]{};
    for (int i = 0; i < 5; ++i) {
        infos[i] = {bufs[i].b, 0, bufs[i].size};
        writes[i].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
        writes[i].dstSet = descriptor;
        writes[i].dstBinding = i;
        writes[i].descriptorCount = 1;
        writes[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        writes[i].pBufferInfo = &infos[i];
    }
    vkUpdateDescriptorSets(device, 5, writes, 0, nullptr);
    VkQueryPoolCreateInfo qp{VK_STRUCTURE_TYPE_QUERY_POOL_CREATE_INFO};
    qp.queryType = VK_QUERY_TYPE_TIMESTAMP;
    qp.queryCount = 4096;
    ck(vkCreateQueryPool(device, &qp, nullptr, &queries));
    VkCommandBufferAllocateInfo ca{VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
    ca.commandPool = pool;
    ca.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    ca.commandBufferCount = 1;
    ck(vkAllocateCommandBuffers(device, &ca, &command));
    VkFenceCreateInfo fi{VK_STRUCTURE_TYPE_FENCE_CREATE_INFO};
    ck(vkCreateFence(device, &fi, nullptr, &fence));
    VkCommandBufferBeginInfo begin{VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    ck(vkBeginCommandBuffer(command, &begin));
    VkBufferCopy rootsCopy{0, 0, bytes[1]};
    vkCmdCopyBuffer(command, rootUpload.b, bufs[1].b, 1, &rootsCopy);
    barrier(command,
            VK_ACCESS_TRANSFER_WRITE_BIT,
            VK_ACCESS_SHADER_READ_BIT,
            VK_PIPELINE_STAGE_TRANSFER_BIT,
            VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT);
    ck(vkEndCommandBuffer(command));
    VkSubmitInfo submit{VK_STRUCTURE_TYPE_SUBMIT_INFO};
    submit.commandBufferCount = 1;
    submit.pCommandBuffers = &command;
    ck(vkQueueSubmit(queue, 1, &submit, fence));
    ck(vkWaitForFences(device, 1, &fence, 1, UINT64_MAX));
    ck(vkResetFences(device, 1, &fence));
    ck(vkResetCommandBuffer(command, 0));
    ck(vkBeginCommandBuffer(command, &begin));
    vkCmdResetQueryPool(command, queries, 0, 2);
    VkBufferCopy copyIn{0, 0, bytes[0]};
    vkCmdCopyBuffer(command, upload.b, bufs[0].b, 1, &copyIn);
    barrier(command,
            VK_ACCESS_TRANSFER_WRITE_BIT,
            VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT,
            VK_PIPELINE_STAGE_TRANSFER_BIT,
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
                        unsigned count = 0) {
        uint32_t push[] = {n, jobs, phase, span, inv[0], inv[1], offset, count};
        vkCmdPushConstants(command, layout, VK_SHADER_STAGE_COMPUTE_BIT, 0, 32, push);
        unsigned groups = (tasks + 255) / 256;
        vkCmdDispatch(command, std::min(65535u, groups), (groups + 65534) / 65535, 1);
        barrier(command,
                VK_ACCESS_SHADER_WRITE_BIT,
                VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT,
                VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
                VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT);
    };
    dispatch(0, jobs * 4ull * n);
    for (unsigned len = n; len > 512; len /= 2)
        dispatch(1, jobs * 2ull * n, len);
    dispatch(5, jobs * 2ull * n);
    dispatch(2, jobs * 2ull * n);
    dispatch(6, jobs * uint64_t(n));
    for (unsigned len = 1024; len <= n; len *= 2)
        dispatch(3, jobs * uint64_t(n), len);
    dispatch(4, jobs * uint64_t(n));
    dispatch(7, jobs * uint64_t(n));
    dispatch(8, jobs * ((blocks + 255) / 256) * 256, 0, n, blocks);
    if (blocks > 256)
        dispatch(8, jobs * 256ull, 0, n + blocks, (blocks + 255) / 256);
    dispatch(9, jobs * uint64_t(n));
    vkCmdWriteTimestamp(command, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, queries, 1);
    barrier(command,
            VK_ACCESS_SHADER_WRITE_BIT,
            VK_ACCESS_TRANSFER_READ_BIT,
            VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
            VK_PIPELINE_STAGE_TRANSFER_BIT);
    VkBufferCopy copyOut{0, 0, bytes[4]};
    vkCmdCopyBuffer(command, bufs[4].b, download.b, 1, &copyOut);
    barrier(command,
            VK_ACCESS_TRANSFER_WRITE_BIT,
            VK_ACCESS_HOST_READ_BIT,
            VK_PIPELINE_STAGE_TRANSFER_BIT,
            VK_PIPELINE_STAGE_HOST_BIT);
    ck(vkEndCommandBuffer(command));
    planSize = n;
    planJobs = jobs;
}
NTTGPU::Timing NTTGPU::multiply(const std::vector<Job> &jobs, const std::function<bool()> &cancel) {
    const auto begin = std::chrono::steady_clock::now();
    if (cancel())
        throw std::runtime_error("NTT cancelled");
    size_t limbs = 1;
    for (auto &j: jobs) {
        if (!j.a || !j.b || !j.result)
            throw std::runtime_error("null integer");
        limbs = std::max(
            {limbs, (mpz_sizeinbase(j.a, 2) + 15) / 16, (mpz_sizeinbase(j.b, 2) + 15) / 16});
    }
    if (limbs > (1u << 21))
        throw std::runtime_error("NTT precision capacity");
    // Each unsigned convolution coefficient <= limbs*(2^16-1)^2.
    if (uint64_t(limbs) * 65535 * 65535 >= uint64_t(moduli[0]) * moduli[1])
        throw std::runtime_error("CRT coefficient bound");
    unsigned n = std::max(512u, std::bit_ceil(unsigned(limbs * 2)));
    plan(n, jobs.size());
    auto *in = static_cast<uint32_t *>(upload.mapped);
    std::memset(in, 0, upload.size);
    std::vector<uint16_t> digits(limbs);
    std::vector<int> signs(jobs.size());
    for (size_t j = 0; j < jobs.size(); ++j) {
        signs[j] = mpz_sgn(jobs[j].a) * mpz_sgn(jobs[j].b);
        for (unsigned op = 0; op < 2; ++op) {
            size_t written = 0;
            mpz_export(digits.data(), &written, -1, 2, 0, 0, op ? jobs[j].b : jobs[j].a);
            for (unsigned p = 0; p < 2; ++p) {
                auto *dst = in + ((j * 2 + p) * 2 + op) * n;
                for (size_t i = 0; i < written; ++i)
                    dst[i] = digits[i];
            }
        }
    }
    if (cancel())
        throw std::runtime_error("NTT cancelled");
    VkSubmitInfo submit{VK_STRUCTURE_TYPE_SUBMIT_INFO};
    submit.commandBufferCount = 1;
    submit.pCommandBuffers = &command;
    ck(vkResetFences(device, 1, &fence));
    auto start = std::chrono::steady_clock::now();
    ck(vkQueueSubmit(queue, 1, &submit, fence));
    bool stopped = false;
    for (;;) {
        auto r = vkWaitForFences(device, 1, &fence, 1, 10000000);
        stopped = stopped || cancel();
        if (r == VK_SUCCESS)
            break;
        if (r != VK_TIMEOUT)
            ck(r);
    }
    Timing timing;
    timing.submit = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    if (stopped)
        throw std::runtime_error("NTT cancelled");
    uint64_t ticks[2];
    ck(vkGetQueryPoolResults(device,
                             queries,
                             0,
                             2,
                             sizeof(ticks),
                             ticks,
                             8,
                             VK_QUERY_RESULT_64_BIT | VK_QUERY_RESULT_WAIT_BIT));
    timing.gpu = (ticks[1] - ticks[0]) * props.limits.timestampPeriod * 1e-9;
    auto *out = static_cast<uint32_t *>(download.mapped);
    // Import only after every GPU task and cancellation check has completed.
    digits.resize(n);
    for (size_t j = 0; j < jobs.size(); ++j) {
        for (size_t i = 0; i < n; ++i)
            digits[i] = uint16_t(out[j * n + i]);
        mpz_import(jobs[j].result, n, -1, 2, 0, 0, digits.data());
        if (signs[j] < 0)
            mpz_neg(jobs[j].result, jobs[j].result);
    }
    timing.total = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
    return timing;
}
