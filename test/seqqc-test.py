# Read generation and an independent FastQC-style oracle for seqqc.
import math
import os
import random
import tempfile

_dir = None
_counts = []
# Sequences followed exactly, as FastQC follows its first 100k.
HEAVY_CAPACITY = 100000
ADAPTERS = ["AGATCGGAAGAG", "CTGTCTCTTATA", "TGGAATTCTCGG", "AAAAAAAAAAAA"]


def tmp_path(name):
    global _dir
    if _dir is None:
        _dir = tempfile.mkdtemp(prefix="bioseqqc-test-")
    return os.path.join(_dir, name)


def make_reads(path, n, seed):
    """n reads of 60-100 bases: 5% copies of one overrepresented sequence,
    20% copies drawn from 50 templates, 10% reading into the Illumina
    adapter, 2% with an N, some lowercase."""
    rng = random.Random(seed)
    over = "".join(rng.choice("ACGT") for _ in range(80))
    templates = ["".join(rng.choice("ACGT") for _ in range(rng.randint(60, 100))) for _ in range(50)]
    with open(path, "w") as f:
        for i in range(n):
            u = rng.random()
            if u < 0.05:
                s = over
            elif u < 0.25:
                s = rng.choice(templates)
            else:
                s = "".join(rng.choice("ACGT") for _ in range(rng.randint(60, 100)))
                if rng.random() < 0.1:
                    k = rng.randint(10, len(s) - 20)
                    s = (s[:k] + "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" + s)[: len(s)]
                if rng.random() < 0.02:
                    k = rng.randrange(len(s))
                    s = s[:k] + "N" + s[k + 1 :]
                if rng.random() < 0.05:
                    s = s.lower()
            q = "".join(chr(33 + rng.randint(2, 41)) for _ in s)
            f.write("@r%d\n%s\n+\n%s\n" % (i, s, q))


def oracle(path):
    with open(path) as f:
        lines = [l.rstrip("\n") for l in f]
    reads = [(lines[i + 1], lines[i + 3]) for i in range(0, len(lines), 4)]
    n = len(reads)
    lengths = [len(s) for s, _ in reads]
    width = max(lengths)
    qs = [[0, 0] for _ in range(width)]
    ns = [[0, 0] for _ in range(width)]
    gcsum, gcn, qtot, btot = 0.0, 0, 0, 0
    starts = [[0] * width for _ in ADAPTERS]
    keys = {}
    for s, q in reads:
        up = s.upper()
        acgt = sum(1 for c in up if c in "ACGT")
        gc = sum(1 for c in up if c in "GC")
        if acgt:
            gcsum += math.floor(100 * gc / acgt + 0.5)
            gcn += 1
        for i, (c, x) in enumerate(zip(up, q)):
            p = min(max(ord(x) - 33, 0), 41)
            qs[i][0] += p
            qs[i][1] += 1
            ns[i][0] += c not in "ACGT"
            ns[i][1] += 1
            qtot += p
            btot += 1
        for a, ad in enumerate(ADAPTERS):
            k = up.find(ad)
            if k >= 0:
                starts[a][k] += 1
        key = up[:50] if len(up) > 75 else up
        keys[key] = keys.get(key, 0) + 1
    content = []
    for a in range(len(ADAPTERS)):
        acc, row = 0, []
        for i in range(width):
            acc += starts[a][i]
            row.append(100.0 * acc / n)
        content.append(row)
    levels = [0.0] * 10
    for c in keys.values():
        levels[min(c, 10) - 1] += 1
    global _counts
    _counts = list(keys.items())
    over = sorted(((k, c, 100.0 * c / n) for k, c in keys.items() if c > 0.001 * n), key=lambda x: (-x[1], x[0]))
    return {
        "reads": n,
        "bases": sum(lengths),
        "minLength": min(lengths),
        "maxLength": max(lengths),
        "meanLength": sum(lengths) / n,
        "gcPercent": gcsum / gcn,
        "meanQuality": qtot / btot,
        "positionQuality": [a / b for a, b in qs],
        "positionN": [100.0 * a / b for a, b in ns],
        "distinct": len(keys),
        "dedupPercent": 100.0 * len(keys) / n,
        "duplicationLevels": [100.0 * x / len(keys) for x in levels],
        "overrepresented": over,
        "adapterContent": content,
    }


def compare(got, want):
    """'' when two reports agree (reals within 1e-9), else the first difference."""
    for k, w in want.items():
        g = got[k]
        if isinstance(w, float):
            if abs(g - w) > 1e-9:
                return "%s: got %r, want %r" % (k, g, w)
        elif isinstance(w, list) and w and isinstance(w[0], list):
            for i, (gr, wr) in enumerate(zip(g, w)):
                if len(gr) != len(wr) or any(abs(x - y) > 1e-9 for x, y in zip(gr, wr)):
                    return "%s[%d]: got %r, want %r" % (k, i, gr[:8], wr[:8])
        elif isinstance(w, list) and w and isinstance(w[0], float):
            if len(g) != len(w) or any(abs(x - y) > 1e-9 for x, y in zip(g, w)):
                return "%s: got %r, want %r" % (k, g[:8], w[:8])
        elif k == "overrepresented" and len(_counts) <= HEAVY_CAPACITY:
            if got["overrepresentedError"] != 0:
                return "overrepresented: error %d with %d distinct sequences" % (got["overrepresentedError"], len(_counts))
            gs = [(s, c) for s, c, _ in g]
            ws = [(s, c) for s, c, _ in w]
            if gs != ws:
                return "overrepresented: got %r, want %r" % (gs[:4], ws[:4])
            if any(abs(x[2] - y[2]) > 1e-9 for x, y in zip(g, w)):
                return "overrepresented: percentages differ"
        elif k == "overrepresented":
            err = got["overrepresentedError"]
            counts = dict(_counts)
            seen = {s for s, _, _ in g}
            for s, c, _ in w:
                if s not in seen:
                    return "overrepresented: %s (%d) missing" % (s, c)
            for s, c, _ in g:
                true = counts.get(s, 0)
                if not (c <= true <= c + err):
                    return "overrepresented: %s reported %d, true %d, error %d" % (s, c, true, err)
                if c + err <= 0.001 * want["reads"]:
                    return "overrepresented: %s cannot reach the threshold" % s
        elif g != w:
            return "%s: got %r, want %r" % (k, g, w)
    return ""


def check(path, report):
    return compare(report, oracle(path))


def fastqc_modules(text):
    """{module name: (flag, header lines, rows)} from fastqc_data.txt text."""
    lines = text.rstrip("\n").split("\n")
    if not lines or not lines[0].startswith("##FastQC\t"):
        raise ValueError("no ##FastQC header")
    mods, cur = {}, None
    for l in lines[1:]:
        if l.startswith(">>END_MODULE"):
            cur = None
        elif l.startswith(">>"):
            name, flag = l[2:].split("\t")
            cur = mods[name] = (flag, [], [])
        elif cur is None:
            raise ValueError("line outside a module: %r" % l)
        elif l.startswith("#"):
            cur[1].append(l)
        else:
            cur[2].append(l.split("\t"))
    return mods


FASTQC_MODULES = [
    "Basic Statistics",
    "Per base sequence quality",
    "Per sequence quality scores",
    "Per base sequence content",
    "Per sequence GC content",
    "Per base N content",
    "Sequence Length Distribution",
    "Sequence Duplication Levels",
    "Overrepresented sequences",
    "Adapter Content",
]


def check_fastqc_data(path, text):
    """'' when fastqc_data.txt text has FastQC's modules and agrees with the
    oracle on the exact columns, else the first difference."""
    try:
        mods = fastqc_modules(text)
    except ValueError as e:
        return str(e)
    if list(mods) != FASTQC_MODULES:
        return "modules: %r" % list(mods)
    for name, (flag, _, _) in mods.items():
        if flag not in ("pass", "warn", "fail"):
            return "%s: flag %r" % (name, flag)
    want = oracle(path)
    basic = dict((r[0], r[1]) for r in mods["Basic Statistics"][2])
    if basic.get("Total Sequences") != str(want["reads"]):
        return "Total Sequences: %r" % basic.get("Total Sequences")
    lens = "%d-%d" % (want["minLength"], want["maxLength"])
    if basic.get("Sequence length") != lens:
        return "Sequence length: %r, want %r" % (basic.get("Sequence length"), lens)
    gc = acgt = 0
    with open(path) as f:
        for i, l in enumerate(f):
            if i % 4 == 1:
                u = l.strip().upper()
                gc += u.count("G") + u.count("C")
                acgt += sum(u.count(c) for c in "ACGT")
    if basic.get("%GC") != str(gc * 100 // acgt):
        return "%%GC: %r, want %d" % (basic.get("%GC"), gc * 100 // acgt)
    rows = mods["Per base sequence quality"][2]
    if len(rows) != want["maxLength"]:
        return "quality rows: %d" % len(rows)
    for i, r in enumerate(rows):
        if int(r[0]) != i + 1 or abs(float(r[1]) - want["positionQuality"][i]) > 1e-5:
            return "quality row %d: %r" % (i, r)
    for i, r in enumerate(mods["Per base N content"][2]):
        if abs(float(r[1]) - want["positionN"][i]) > 1e-5:
            return "N row %d: %r" % (i, r)
    counts = {}
    with open(path) as f:
        for i, l in enumerate(f):
            if i % 4 == 1:
                n = len(l.rstrip("\n"))
                counts[n] = counts.get(n, 0) + 1
    got = dict((int(r[0]), float(r[1])) for r in mods["Sequence Length Distribution"][2])
    if got != {k: float(v) for k, v in counts.items()}:
        return "lengths differ"
    dedup = float(mods["Sequence Duplication Levels"][1][0].split("\t")[1])
    if abs(dedup - want["dedupPercent"]) > 1e-5:
        return "dedup: %r, want %r" % (dedup, want["dedupPercent"])
    over = mods["Overrepresented sequences"][2]
    if [r[0] for r in over][:1] != [s for s, _, _ in want["overrepresented"]][:1]:
        return "top overrepresented: %r" % over[:1]
    if len(mods["Adapter Content"][2]) == 0:
        return "no adapter rows"
    return ""
