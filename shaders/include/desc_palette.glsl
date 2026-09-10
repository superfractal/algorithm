// Modified by GPT-6 on 2026-09-07
#ifndef DESC_PALETTE_INCLUDE
#define DESC_PALETTE_INCLUDE

layout (set = DESC_PALETTE, binding = 0) buffer PaletteSSBO {
    uint size;
    float interval;
    double offset;
    uint coloring;
    uint single_coloring;
    float animation_speed;
    uint linear_interpolation;
    vec4 palette[];
} palette_settings;

#endif
