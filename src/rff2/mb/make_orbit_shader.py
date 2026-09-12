# Created by GPT-6 on 2026-09-11
# Modified by GPT-6 on 2026-09-12
"""Generate the eight orbit shaders from readable GLSL source fragments.

Replacement anchors are GLSL token sequences, so changing indentation does not
disable an optimization. Each replacement checks its expected occurrence count
and fails before writing output if the source structure no longer matches.
"""

from pathlib import Path
import re
import textwrap


# Comments and whitespace are ignored for matching; every other character is
# retained in a token. Multi-character operators must remain indivisible.
GLSL_TOKEN = re.compile(
    r"(?P<comment>//[^\n]*|/\*[\s\S]*?\*/)"
    r"|[A-Za-z_]\w*"
    r"|(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?[uUlLfF]*"
    r"|<<=|>>=|\+\+|--|==|!=|<=|>=|&&|\|\||<<|>>"
    r"|[+*/%&|^!-]=|[^\s]"
)


def token_spans(source):
    """Return non-comment tokens together with their original text positions."""
    return [
        match for match in GLSL_TOKEN.finditer(source)
        if match.lastgroup != "comment"
    ]


def anchor_spans(source, anchor, expected=1):
    tokens = token_spans(source)
    values = [token.group() for token in tokens]
    wanted = [token.group() for token in token_spans(anchor)]
    if not wanted:
        raise ValueError("A shader transformation anchor must not be empty")

    spans = [
        (tokens[index].start(), tokens[index + len(wanted) - 1].end())
        for index in range(len(tokens) - len(wanted) + 1)
        if values[index:index + len(wanted)] == wanted
    ]
    if len(spans) != expected:
        raise ValueError(
            f"Shader anchor expected {expected} occurrence(s), found {len(spans)}: "
            f"{anchor}"
        )
    return spans


def replace_tokens(source, anchor, replacement, expected=1):
    """Replace from right to left so earlier source positions stay valid."""
    for start, end in reversed(anchor_spans(source, anchor, expected)):
        if "\n" in replacement:
            line_start = source.rfind("\n", 0, start) + 1
            indentation = source[line_start:start]
            if not indentation.isspace():
                indentation = ""
            lines = textwrap.dedent(replacement).strip().splitlines()
            rendered = lines[0] + "".join(
                "\n" + (line if line.startswith("#") else indentation + line)
                for line in lines[1:]
            )
        else:
            rendered = replacement
        source = source[:start] + rendered + source[end:]
    return source


def indent_glsl(source, prefix):
    """Indent phase code while leaving preprocessor directives at column one."""
    return "\n".join(
        line if line.startswith("#") else prefix + line
        for line in source.splitlines()
    )


COMPLEX_PRODUCT_PHASE = """\
if (p.phase == 2) {
    if (t >= (p.jobs / PRODUCTS) * 2 * p.n)
        return;
    uint stream = t / p.n, seg = stream / 2, prime = stream % 2, i = t % p.n,
         base = seg * PRODUCTS * 4;
    uint x = data[(base + prime * 2) * p.n + i], y = data[(base + 4 + prime * 2) * p.n + i];
#if COMPLEX_SQUARE
    uint st = p.jobs * p.n, L = p.n / 2, meta = st + 3 * (p.jobs / PRODUCTS) * p.n,
         at = result[meta + 12 * seg + 4];
    bool negative =
        at != 0 && result[st + 2 * seg * L + at - 1] < result[st + (2 * seg + 1) * L + at - 1];
    uint xx = mont(x, x, prime), yy = mont(y, y, prime), xy = mont(x, y, prime);
    data[(base + prime * 2) * p.n + i] = negative ? minus(yy, xx, prime) : minus(xx, yy, prime);
    data[(base + 4 + prime * 2) * p.n + i] = plus(xy, xy, prime);
#else
    data[(base + prime * 2) * p.n + i] = mont(x, x, prime);
    data[(base + 4 + prime * 2) * p.n + i] = mont(y, y, prime);
    data[(base + 8 + prime * 2) * p.n + i] = mont(x, y, prime);
#endif
    return;
}
"""


SIGNED_CRT_RECONSTRUCTION = """\
uint64_t value = uint64_t(a) + uint64_t(2013265921u) * uint64_t(d);
#if COMPLEX_SQUARE
const uint64_t Q = uint64_t(2013265921u) * uint64_t(1811939329u);
if (value > Q / 2)
    value -= Q;
#endif
coefficients[t] = value;
return;
"""


ORBIT_HELPERS = """\
shared uvec4 maxima[256];
uvec4 blockMax(uvec4 value, uint lane) {
    maxima[lane] = value;
    barrier();
    for (uint distance = 128; distance > 0; distance /= 2) {
        if (lane < distance)
            maxima[lane] = max(maxima[lane], maxima[lane + distance]);
        barrier();
    }
    return maxima[0];
}
uint signedCompose(uint after, uint before) {
    uint v = 0;
    for (uint c = 0; c < 5; ++c)
        v |= ((after >> (((before >> (3 * c)) & 7) * 3)) & 7) << (3 * c);
    return v;
}
uint signedCarry(uint m, uint c) {
    return (m >> (3 * c)) & 7;
}
void signedScan(uint lane) {
    for (uint step = 1; step < 256; step *= 2) {
        uint prev = lane >= step ? tile[lane - step] : 18056u;
        barrier();
        tile[lane] = signedCompose(tile[lane], prev);
        barrier();
    }
}
"""


ORIGINAL_MAXIMUM_UPDATES = """\
if (x != y)
    atomicMax(result[meta + 12 * seg + 4], i + 1);
if (x != 0)
    atomicMax(result[meta + 12 * seg + 8], i + 1);
if (y != 0)
    atomicMax(result[meta + 12 * seg + 9], i + 1);
"""


GROUPED_MAXIMUM_UPDATES = """\
#if GROUP_REDUCE
uvec4 highest =
    blockMax(uvec4(x != y ? i + 1 : 0, x != 0 ? i + 1 : 0, y != 0 ? i + 1 : 0, 0), lane);
if (lane == 0) {
    atomicMax(result[meta + 12 * seg + 4], highest.x);
    atomicMax(result[meta + 12 * seg + 8], highest.y);
    atomicMax(result[meta + 12 * seg + 9], highest.z);
}
#else
if (x != y)
    atomicMax(result[meta + 12 * seg + 4], i + 1);
if (x != 0)
    atomicMax(result[meta + 12 * seg + 8], i + 1);
if (y != 0)
    atomicMax(result[meta + 12 * seg + 9], i + 1);
#endif
"""


GROUPED_OVERFLOW_UPDATES = """\
#if GROUP_REDUCE
uvec4 highest = blockMax(uvec4(v != result[cs + t] ? i + 1 : 0, lost != 0 ? 1 : 0, 0, 0), lane);
if (lane == 0) {
    atomicMax(result[meta + 12 * seg + 5 + c], highest.x);
    atomicOr(result[meta + 12 * seg + 7], highest.y);
}
#else
if (v != result[cs + t])
    atomicMax(result[meta + 12 * seg + 5 + c], i + 1);
if (lost != 0)
    atomicOr(result[meta + 12 * seg + 7], 1u);
#endif
return;
"""


def generate_shader(source_directory):
    """Construct the shared shader body without changing arithmetic order."""
    source = (source_directory / "ntt.comp").read_text()

    # Orbit batches have two inputs per segment but two or three products.
    source = replace_tokens(
        source,
        "uint streams=p.jobs*(inverse?2:4);",
        "uint streams = inverse ? p.jobs * 2 : (p.jobs / PRODUCTS) * 4;",
    )
    source = replace_tokens(
        source,
        "uint streams=p.jobs*(inverse?2:4),blocks=p.n/512;",
        "uint streams = inverse ? p.jobs * 2 : (p.jobs / PRODUCTS) * 4, "
        "blocks = p.n / 512;",
    )
    source = replace_tokens(
        source,
        "base=(inverse?stream*2:stream)*p.n",
        "base =\n"
        "            (inverse ? stream * 2 :\n"
        "             ((stream / 4) * PRODUCTS * 4 +\n"
        "              (stream % 4) / 2 * 2 + (stream % 2) * 4)) * p.n",
        expected=2,
    )
    source = replace_tokens(
        source,
        "prime=inverse?stream%2:(stream/2)%2",
        "prime = inverse ? stream % 2 : (stream % 4) / 2",
        expected=2,
    )

    # Replace the ordinary product and inverse-transform phase together.
    start = anchor_spans(source, "if (p.phase == 2) {")[0][0]
    end = anchor_spans(source, "if (p.phase == 4) {")[0][0]
    if end <= start:
        raise ValueError("Shader product phase anchors are out of order")
    product_phase = indent_glsl(COMPLEX_PRODUCT_PHASE.strip(), "    ")
    source = source[:start] + product_phase.lstrip() + "\n    " + source[end:]
    source = replace_tokens(
        source,
        "coefficients[t]=uint64_t(a)+uint64_t(2013265921u)*uint64_t(d);return;",
        SIGNED_CRT_RECONSTRUCTION,
    )
    source = replace_tokens(source, "void main(){", ORBIT_HELPERS + "\nvoid main() {")

    phases = (source_directory / "orbit_phases.inc").read_text()
    phases = replace_tokens(phases, "p.jobs/3", "p.jobs / PRODUCTS")
    phases = replace_tokens(phases, "seg*12", "seg * PRODUCTS * 4", expected=2)
    phases = replace_tokens(phases, "3*seg", "PRODUCTS * seg", expected=6)
    phases = replace_tokens(
        phases,
        "uint shift=p.offset-c,source=(PRODUCTS*seg+(c==0?0:2))*p.n",
        "uint shift = p.offset - (COMPLEX_SQUARE == 1 ? 0 : c),\n"
        "     source =\n"
        "         (PRODUCTS * seg + (c == 0 ? 0 : (COMPLEX_SQUARE == 1 ? 1 : 2))) * p.n",
    )
    phases = replace_tokens(phases, ORIGINAL_MAXIMUM_UPDATES, GROUPED_MAXIMUM_UPDATES)
    phases = replace_tokens(
        phases,
        "result[tmp+t]=v;if(v!=result[cs+t])atomicMax(result[meta+12*seg+5+c],i+1);",
        "result[tmp + t] = v;",
    )
    phases = replace_tokens(
        phases,
        "if(lost!=0)atomicOr(result[meta+12*seg+7],1u);return;",
        GROUPED_OVERFLOW_UPDATES,
    )

    # Both include fragments contain main() phase branches, not declarations.
    closing = token_spans(source)[-1]
    if closing.group() != "}":
        raise ValueError("Shader main() closing brace is missing")
    optimized = (source_directory / "optimized_phases.inc").read_text()
    appended = indent_glsl(optimized.rstrip() + "\n" + phases.rstrip(), "    ")
    return source[:closing.start()] + appended + "\n}\n"


def write_variants(source_directory):
    source = generate_shader(source_directory)
    variants = {}
    for mask in range(8):
        definitions = (
            "#version 460\n"
            f"#define GROUP_REDUCE {int(bool(mask & 1))}\n"
            f"#define COMPLEX_SQUARE {int(bool(mask & 2))}\n"
            f"#define PRODUCTS {2 if mask & 2 else 3}"
        )
        header = (
            "// Created by GPT-6 on 2026-09-12\n"
            "// Modified by GPT-6 on 2026-09-12\n"
            "// Generated variant; original source notices follow.\n"
        )
        variants[mask] = header + replace_tokens(source, "#version 460", definitions)

    # Validate every transformation before replacing any generated file.
    for mask, variant in variants.items():
        (source_directory / f"orbit_{mask}.comp").write_text(variant)


if __name__ == "__main__":
    write_variants(Path(__file__).resolve().parent)
    print("Generated eight shader variants")
