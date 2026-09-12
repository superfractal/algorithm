// Created by GPT-6 on 2026-09-11
// Vulkan NTT runtime, adapted from the local v197 device/buffer setup.
// Modified by GPT-6 on 2026-09-12
#pragma once
#include <functional>
#include <limits>
#include <vulkan/vulkan.h>
#include <gmp.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>
#include <cstring>
#include <stdexcept>
#include <algorithm>
#include <omp.h>
inline void ck(VkResult r) {
    if (r != VK_SUCCESS)
        throw std::runtime_error("Vulkan error " + std::to_string(r));
}
struct Z {
    mpz_t x;
    Z() {
        mpz_init(x);
    }
    Z(const Z &o) {
        mpz_init_set(x, o.x);
    }
    ~Z() {
        mpz_clear(x);
    }
};
struct Buffer {
    VkBuffer b{};
    VkDeviceMemory m{};
    void *mapped{};
    VkDeviceSize size;
};
struct NTTGPU {
    VkInstance instance{};
    VkPhysicalDevice physical{};
    VkDevice device{};
    VkQueue queue{};
    uint32_t family = 0;
    VkDescriptorSetLayout setLayout{};
    VkPipelineLayout layout{};
    VkPipeline pipeline{};
    VkDescriptorPool descriptors{};
    VkCommandPool pool{};
    VkPhysicalDeviceProperties props{};
    NTTGPU(const char *shader) {
        VkApplicationInfo app{VK_STRUCTURE_TYPE_APPLICATION_INFO};
        app.apiVersion = VK_API_VERSION_1_2;
        VkInstanceCreateInfo ci{VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
        ci.pApplicationInfo = &app;
        ck(vkCreateInstance(&ci, nullptr, &instance));
        uint32_t n = 0;
        ck(vkEnumeratePhysicalDevices(instance, &n, nullptr));
        std::vector<VkPhysicalDevice> devices(n);
        ck(vkEnumeratePhysicalDevices(instance, &n, devices.data()));
        for (auto d: devices) {
            VkPhysicalDeviceProperties p;
            vkGetPhysicalDeviceProperties(d, &p);
            if (!physical || p.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU)
                physical = d;
        }
        if (!physical)
            throw std::runtime_error("no Vulkan GPU");
        vkGetPhysicalDeviceProperties(physical, &props);
        VkPhysicalDeviceFeatures features;
        vkGetPhysicalDeviceFeatures(physical, &features);
        if (!features.shaderInt64)
            throw std::runtime_error("shaderInt64 unavailable");
        uint32_t qn;
        vkGetPhysicalDeviceQueueFamilyProperties(physical, &qn, nullptr);
        std::vector<VkQueueFamilyProperties> qs(qn);
        vkGetPhysicalDeviceQueueFamilyProperties(physical, &qn, qs.data());
        while (family < qn && !(qs[family].queueFlags & VK_QUEUE_COMPUTE_BIT))
            ++family;
        if (family == qn || !qs[family].timestampValidBits)
            throw std::runtime_error("compute timestamp queue unavailable");
        float priority = 1;
        VkDeviceQueueCreateInfo qi{VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
        qi.queueFamilyIndex = family;
        qi.queueCount = 1;
        qi.pQueuePriorities = &priority;
        VkPhysicalDeviceFeatures enabled{};
        enabled.shaderInt64 = 1;
        VkDeviceCreateInfo di{VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
        di.queueCreateInfoCount = 1;
        di.pQueueCreateInfos = &qi;
        di.pEnabledFeatures = &enabled;
        ck(vkCreateDevice(physical, &di, nullptr, &device));
        vkGetDeviceQueue(device, family, 0, &queue);
        VkDescriptorSetLayoutBinding bindings[5]{};
        for (int i = 0; i < 5; ++i) {
            bindings[i].binding = i;
            bindings[i].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
            bindings[i].descriptorCount = 1;
            bindings[i].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
        }
        VkDescriptorSetLayoutCreateInfo sl{VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
        sl.bindingCount = 5;
        sl.pBindings = bindings;
        ck(vkCreateDescriptorSetLayout(device, &sl, nullptr, &setLayout));
        VkPushConstantRange push{VK_SHADER_STAGE_COMPUTE_BIT, 0, 32};
        VkPipelineLayoutCreateInfo pl{VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
        pl.setLayoutCount = 1;
        pl.pSetLayouts = &setLayout;
        pl.pushConstantRangeCount = 1;
        pl.pPushConstantRanges = &push;
        ck(vkCreatePipelineLayout(device, &pl, nullptr, &layout));
        std::ifstream in(shader, std::ios::binary | std::ios::ate);
        if (!in)
            throw std::runtime_error("missing SPIR-V");
        size_t bytes = in.tellg();
        std::vector<uint32_t> code((bytes + 3) / 4);
        in.seekg(0);
        in.read((char *) code.data(), bytes);
        VkShaderModuleCreateInfo sm{VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
        sm.codeSize = bytes;
        sm.pCode = code.data();
        VkShaderModule module;
        ck(vkCreateShaderModule(device, &sm, nullptr, &module));
        VkComputePipelineCreateInfo pc{VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO};
        pc.stage.sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
        pc.stage.stage = VK_SHADER_STAGE_COMPUTE_BIT;
        pc.stage.module = module;
        pc.stage.pName = "main";
        pc.layout = layout;
        ck(vkCreateComputePipelines(device, {}, 1, &pc, nullptr, &pipeline));
        vkDestroyShaderModule(device, module, nullptr);
        VkDescriptorPoolSize ps{VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 5};
        VkDescriptorPoolCreateInfo dp{VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
        dp.maxSets = 1;
        dp.poolSizeCount = 1;
        dp.pPoolSizes = &ps;
        ck(vkCreateDescriptorPool(device, &dp, nullptr, &descriptors));
        VkCommandPoolCreateInfo cp{VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO};
        cp.queueFamilyIndex = family;
        cp.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
        ck(vkCreateCommandPool(device, &cp, nullptr, &pool));
        std::cout << "GPU=" << props.deviceName
                  << " timestampPeriod=" << props.limits.timestampPeriod << std::endl;
    }
    Buffer buffer(VkDeviceSize size, bool host) {
        Buffer b{};
        b.size = size;
        VkBufferCreateInfo bi{VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
        bi.size = size;
        bi.usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_TRANSFER_SRC_BIT |
                   VK_BUFFER_USAGE_TRANSFER_DST_BIT;
        ck(vkCreateBuffer(device, &bi, nullptr, &b.b));
        VkMemoryRequirements req;
        vkGetBufferMemoryRequirements(device, b.b, &req);
        VkPhysicalDeviceMemoryProperties mp;
        vkGetPhysicalDeviceMemoryProperties(physical, &mp);
        const VkMemoryPropertyFlags flags =
            host ? VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT
                 : VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT;
        uint32_t t = 0;
        while (t < mp.memoryTypeCount && (!(req.memoryTypeBits & (1u << t)) ||
                                          (mp.memoryTypes[t].propertyFlags & flags) != flags))
            ++t;
        if (t == mp.memoryTypeCount)
            throw std::runtime_error("memory type unavailable");
        VkMemoryAllocateInfo ai{VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
        ai.allocationSize = req.size;
        ai.memoryTypeIndex = t;
        ck(vkAllocateMemory(device, &ai, nullptr, &b.m));
        ck(vkBindBufferMemory(device, b.b, b.m, 0));
        if (host)
            ck(vkMapMemory(device, b.m, 0, size, 0, &b.mapped));
        return b;
    }
    void free(Buffer &b) {
        if (b.mapped)
            vkUnmapMemory(device, b.m);
        vkDestroyBuffer(device, b.b, nullptr);
        vkFreeMemory(device, b.m, nullptr);
    }
    void barrier(VkCommandBuffer cmd,
                 VkAccessFlags src,
                 VkAccessFlags dst,
                 VkPipelineStageFlags from,
                 VkPipelineStageFlags to) {
        VkMemoryBarrier b{VK_STRUCTURE_TYPE_MEMORY_BARRIER};
        b.srcAccessMask = src;
        b.dstAccessMask = dst;
        vkCmdPipelineBarrier(cmd, from, to, 0, 1, &b, 0, nullptr, 0, nullptr);
    }

    static uint32_t power(uint32_t a, uint64_t n, uint32_t p) {
        uint64_t v = 1;
        for (; n; n >>= 1, a = uint64_t(a) * a % p)
            if (n & 1)
                v = v * a % p;
        return uint32_t(v);
    }
    struct Timing {
        double gpu = 0, submit = 0, total = 0;
    };
    struct Job {
        mpz_srcptr a, b;
        mpz_ptr result;
    };
    static constexpr uint32_t moduli[2] = {2013265921, 1811939329};
    uint32_t size = 0;
    std::vector<uint32_t> rootTable;
    void prepare(unsigned n) {
        if (n == size)
            return;
        if (n < 512 || n > (1u << 22) || (n & (n - 1)))
            throw std::runtime_error("NTT length unsupported");
        rootTable.resize(2ull * n);
        const uint32_t generators[2] = {31, 13};
        for (unsigned p = 0; p < 2; ++p) {
            const auto q = moduli[p];
            if ((q - 1) % n)
                throw std::runtime_error("NTT modulus capacity");
            const auto w = power(generators[p], (q - 1) / n, q), inv = power(w, q - 2, q);
            if (power(w, n, q) != 1 || power(w, n / 2, q) == 1)
                throw std::runtime_error("NTT root order");
            for (unsigned inverse = 0; inverse < 2; ++inverse) {
                uint64_t v = (1ull << 32) % q;
                for (unsigned j = 0; j < n / 2; ++j) {
                    rootTable[p * n + inverse * n / 2 + j] = v;
                    v = v * (inverse ? inv : w) % q;
                }
            }
        }
        size = n;
    }
    Buffer bufs[5]{}, upload{}, download{}, rootUpload{};
    VkDescriptorSet descriptor{};
    VkCommandBuffer command{};
    VkQueryPool queries{};
    VkFence fence{};
    unsigned planSize = 0, planJobs = 0;
    void resetPlan();
    void plan(unsigned n, unsigned jobs);
    Timing multiply(
        const std::vector<Job> &jobs, const std::function<bool()> &cancel = [] {
            return false;
        });
    ~NTTGPU() {
        if (device) {
            vkDeviceWaitIdle(device);
            resetPlan();
            vkDestroyCommandPool(device, pool, nullptr);
            vkDestroyDescriptorPool(device, descriptors, nullptr);
            vkDestroyPipeline(device, pipeline, nullptr);
            vkDestroyPipelineLayout(device, layout, nullptr);
            vkDestroyDescriptorSetLayout(device, setLayout, nullptr);
            vkDestroyDevice(device, nullptr);
        }
        if (instance)
            vkDestroyInstance(instance, nullptr);
    }
};
