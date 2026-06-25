#!/usr/bin/env python3
"""Deterministic synthetic ChIP-seq fixture for the MACS3 model-building stage.

200 peaks 5 kb apart on a 1 Mb `chr1`; each peak has 30 plus-strand reads whose
5' ends cluster ~frag/2 upstream of the summit and 30 minus-strand reads whose
5' ends cluster ~frag/2 downstream, giving the bimodal +/- pattern PeakModel
detects. Plus 4000 uniform background reads. Seeded for reproducibility.

    python3 simulate_chip.py > sim.sam
    samtools sort -o sim.bam sim.sam && samtools index sim.bam

macs3 3.0.4 `predictd -g 1000000` on the result: 125 paired peaks, predicted
fragment length d = 145 (alternatives 29,86,145,204,275). The Rust PeakModel
must reproduce d = 145.
"""
import random
import sys

random.seed(20260625)
L, NPEAKS, FRAG, RL, RPS, BG = 1_000_000, 200, 147, 50, 30, 4000
SEQ, QUAL = "A" * RL, "I" * RL
out, qn = [], 0


def emit(flag, leftmost):
    global qn
    if 0 <= leftmost <= L - RL:
        qn += 1
        out.append(
            f"r{qn}\t{flag}\tchr1\t{leftmost + 1}\t30\t{RL}M\t*\t0\t0\t{SEQ}\t{QUAL}"
        )


for i in range(NPEAKS):
    s = 2500 + i * 5000
    # Distinct 5' positions per read (offset -RPS//2..RPS//2) so that --keep-dup 1
    # does not collapse the signal — callpeak builds its model on deduped tags.
    for j in range(RPS):
        off = j - RPS // 2
        emit(0, s - FRAG // 2 + off)
        emit(16, s + FRAG // 2 + off - RL)
for _ in range(BG):
    emit(random.choice([0, 16]), random.randint(0, L - RL))

sys.stdout.write("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000000\n")
sys.stdout.write("\n".join(out) + "\n")
