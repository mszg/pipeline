"""Conservative read-level support for normalized small variants.

Indels require aligned sequence spanning both the VCF anchor and a right flank.
The window includes equivalent right-shifted placements in a tandem repeat, so
an alternate event shifted in the CIGAR is not mistaken for reference support.
Incomplete/ambiguous/conflicting windows are OTHER, never exact REF or ALT.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class AlleleWindow:
    start: int
    end: int
    ref: str
    alt: str


def indel_window(sequence, start, ref, alt):
    """Return a 0-based half-open, repeat-aware window, or None if unsupported.

    Input alleles must be normalized, left-anchored pure insertions/deletions.
    Both window endpoints must be observed as aligned read bases. Variants at
    a reference end without a right flank cannot be assessed this way.
    """
    sequence, ref, alt = sequence.upper(), ref.upper(), alt.upper()
    if start < 0 or sequence[start:start + len(ref)] != ref:
        raise ValueError("Variant REF does not match the reference sequence")
    if not ref or not alt or set(ref + alt) - set("ACGT"):
        return None
    if len(ref) == 1 and len(alt) > 1 and alt.startswith(ref):
        motif = alt[1:]
    elif len(alt) == 1 and len(ref) > 1 and ref.startswith(alt):
        motif = ref[1:]
    else:
        return None
    end = start + len(ref)
    while end < len(sequence) and sequence[end] == motif[0]:
        motif = motif[1:] + sequence[end]
        end += 1
    if end >= len(sequence):
        return None
    end += 1  # aligned right flank beyond every equivalent repeat placement
    ref_window = sequence[start:end]
    alt_window = alt + sequence[start + len(ref):end]
    if set(ref_window + alt_window) - set("ACGT"):
        return None
    return AlleleWindow(start, end, ref_window, alt_window)


def aligned_window(read, start, end):
    """Extract query sequence inside reference boundaries using the CIGAR.

    Return None unless both boundary bases align (M/=/X) and no reference skip,
    internal clip or padding intersects the window. Insertions inside the
    boundaries are included; insertions immediately outside are excluded.
    """
    if not read.query_sequence or not read.cigartuples or read.is_unmapped:
        return None
    if read.reference_start > start or read.reference_end is None or read.reference_end < end:
        return None
    sequence = read.query_sequence
    rpos, qpos = read.reference_start, 0
    left_aligned = right_aligned = False
    pieces = []
    for op, length in read.cigartuples:
        if op in (0, 7, 8):  # M, =, X
            left_aligned |= rpos <= start < rpos + length
            right_aligned |= rpos <= end - 1 < rpos + length
            lo, hi = max(start, rpos), min(end, rpos + length)
            if lo < hi:
                pieces.append(sequence[qpos + lo - rpos:qpos + hi - rpos])
            rpos += length
            qpos += length
        elif op == 1:  # I
            if start < rpos < end:
                pieces.append(sequence[qpos:qpos + length])
            qpos += length
        elif op in (2, 3):  # D, N
            if op == 3 and rpos < end and rpos + length > start:
                return None
            rpos += length
        elif op in (4, 5, 6):  # S, H, P
            if start < rpos < end:
                return None
            if op == 4:
                qpos += length
        else:
            return None
        if rpos >= end:
            break
    if not (left_aligned and right_aligned):
        return None
    return "".join(pieces).upper()


def classify_indel(read, window):
    observed = aligned_window(read, window.start, window.end)
    if observed == window.ref:
        return "REF"
    if observed == window.alt:
        return "ALT"
    return "OTHER"


def count_indel_support(bam, chrom, window):
    """Count all eligible reads overlapping the canonical anchor, without a cap.

    Mirror mpileup's QC-fail/duplicate exclusions, and require primary mapped
    alignments. BAMs have already passed workflow MAPQ and length filtering.
    """
    counts = {"REF": 0, "ALT": 0, "OTHER": 0}
    for read in bam.fetch(chrom, window.start, window.start + 1):
        if read.flag & (4 | 256 | 512 | 1024 | 2048):
            continue
        counts[classify_indel(read, window)] += 1
    return counts["REF"], counts["ALT"], counts["OTHER"], True


def observations(bases):
    """Parse samtools pileup observations; malformed input fails explicitly."""
    out, i = [], 0
    while i < len(bases):
        c = bases[i]
        if c == "^":
            if i + 1 >= len(bases):
                raise ValueError("Truncated pileup read-start marker")
            i += 2
            continue
        if c == "$":
            i += 1
            continue
        if c not in ".,ACGTNacgtn*#<>":
            raise ValueError(f"Unexpected pileup character {c!r}")
        obs = {"base": c, "indels": []}
        i += 1
        while i < len(bases) and bases[i] in "+-":
            sign = bases[i]
            i += 1
            j = i
            while j < len(bases) and bases[j].isdigit():
                j += 1
            if j == i:
                raise ValueError("Pileup indel is missing its length")
            length = int(bases[i:j])
            if length < 1 or j + length > len(bases):
                raise ValueError("Invalid/truncated pileup indel")
            allele = bases[j:j + length].upper()
            if set(allele) - set("ACGTN*#"):
                raise ValueError("Invalid pileup indel sequence")
            obs["indels"].append((sign, length, allele))
            i = j + length
        out.append(obs)
    return out


def count_snv_support(ref, alt, bases):
    ref, alt = ref.upper(), alt.upper()
    if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT" or ref == alt:
        return 0, 0, len(observations(bases)), False
    ref_n = alt_n = other = 0
    for obs in observations(bases):
        base = obs["base"].upper()
        if obs["indels"]:
            other += 1
        elif base in ".," or base == ref:
            ref_n += 1
        elif base == alt:
            alt_n += 1
        else:
            other += 1
    return ref_n, alt_n, other, True
