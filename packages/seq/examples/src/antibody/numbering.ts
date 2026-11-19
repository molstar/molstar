/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The schemes' rules, in code (examples.md, Antibody): a label a position of a chain, `null` outside the variable
// domain. IMGT numbers by fixed positions: each CDR's positions filled from its ends inward, the gaps at the top of the
// loop; a CDR3 longer than 13 gains `111.1`, `112.1`, `111.2`, … about its middle, so `112.1` comes before `112` (EC-7).
// Kabat is derived from IMGT by counting between the domain's anchors (the cysteines at IMGT 23 and 104, the
// tryptophan at 41 and the phenylalanine or tryptophan at 118, which Kabat numbers L23, L35, L88, L98 and H22, H36,
// H92, H103): each stretch between two anchors is its Kabat labels, a CDR's insertions lettered after its insertion
// point (L27, L95; H35, H52, H100) and its deletions taken from that point down. The rules are checked against the
// wild type's published IMGT numbering and the Kabat labels the literature gives trastuzumab (52A, 82A–C, 100A–C).

export type Chain = 'light' | 'heavy';
export type Labels = readonly (string | null)[];

/** A CDR by IMGT's definition: its first and last position. */
export const ImgtCdrs = { 1: [27, 38], 2: [56, 65], 3: [105, 117] } as const;

/**
 * The order IMGT leaves a shorter CDR's positions empty in (the top of the loop first, the C-terminal side of its
 * middle before the N-terminal one, so an odd count keeps one more position at the start), and the order a CDR3 longer
 * than 13 gains positions in.
 */
const ImgtGaps: Readonly<Record<1 | 2 | 3, readonly number[]>> = {
  1: [33, 32, 34, 31, 35, 30, 36, 29, 37, 28],
  2: [61, 60, 62, 59, 63, 58, 64, 57],
  3: [111, 112, 110, 113, 109, 114, 108, 115, 107, 116],
};

/** IMGT's labels of a CDR of `length` residues, in sequence order. */
export function imgtCdrLabels(cdr: 1 | 2 | 3, length: number): string[] {
  const [first, last] = ImgtCdrs[cdr];
  const size = last - first + 1;
  if (length <= size) {
    const gaps = new Set(ImgtGaps[cdr].slice(0, size - length));
    if (gaps.size !== size - length) throw new Error(`IMGT CDR${cdr} of ${length} residues: too short`);
    const out: string[] = [];
    for (let p = first; p <= last; p++) if (!gaps.has(p)) out.push(`${p}`);
    return out;
  }
  if (cdr !== 3) throw new Error(`IMGT CDR${cdr} of ${length} residues: longer than its ${size} positions`);
  // 111.1 first, then 112.1, 111.2, 112.2, …: 111's run rises, 112's falls, so 112.n comes before 112
  const extra = length - size;
  const left = Math.ceil(extra / 2);
  const right = extra - left;
  return [
    ...Array.from({ length: 111 - first + 1 }, (_, k) => `${first + k}`),
    ...Array.from({ length: left }, (_, k) => `111.${k + 1}`),
    ...Array.from({ length: right }, (_, k) => `112.${right - k}`),
    ...Array.from({ length: last - 112 + 1 }, (_, k) => `${112 + k}`),
  ];
}

/** An IMGT label's place in sequence order: its integer, then 111's insertions rising and 112's falling before 112. */
export function imgtOrder(label: string): number {
  const [int, ins] = label.split('.');
  const n = Number(int);
  if (ins === undefined) return n;
  return n === 112 ? n - Number(ins) / 1000 : n + Number(ins) / 1000;
}

/** A stretch of Kabat labels between two anchors: fixed labels, then a CDR's run with its insertion point, then fixed. */
interface Stretch {
  readonly before: readonly string[];
  readonly from: number;
  readonly to: number;
  /** Insertions are lettered after it; deletions are taken from it down. */
  readonly at: number;
  readonly after: readonly string[];
}

const run = (from: number, to: number) => Array.from({ length: to - from + 1 }, (_, k) => `${from + k}`);

/** Kabat's anchors (their IMGT positions and Kabat labels) and the stretches between them, by chain. */
const Kabat: Readonly<Record<Chain, { readonly anchors: readonly [number, number][]; readonly stretches: Stretch[] }>> =
  {
    light: {
      anchors: [
        [23, 23],
        [41, 35],
        [104, 88],
        [118, 98],
      ],
      stretches: [
        { before: [], from: 24, to: 34, at: 27, after: [] },
        // CDR-L2 is seven residues in both chains the dataset holds; its insertion point is never met
        { before: run(36, 49), from: 50, to: 56, at: 54, after: run(57, 87) },
        { before: [], from: 89, to: 97, at: 95, after: [] },
      ],
    },
    heavy: {
      anchors: [
        [23, 22],
        [41, 36],
        [104, 92],
        [118, 103],
      ],
      stretches: [
        { before: run(23, 30), from: 31, to: 35, at: 35, after: [] },
        {
          before: run(37, 49),
          from: 50,
          to: 65,
          at: 52,
          after: [...run(66, 82), '82A', '82B', '82C', ...run(83, 91)],
        },
        { before: ['93', '94'], from: 95, to: 102, at: 100, after: [] },
      ],
    },
  };

const Letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

/** A CDR's Kabat labels for `length` residues: insertions lettered after `at`, deletions from `at` down. */
function kabatRun(s: Stretch, length: number): string[] {
  const base = run(s.from, s.to);
  if (length >= base.length) {
    const k = base.indexOf(`${s.at}`) + 1;
    const ins = Array.from({ length: length - base.length }, (_, i) => `${s.at}${Letters[i]}`);
    return [...base.slice(0, k), ...ins, ...base.slice(k)];
  }
  const drop = new Set(Array.from({ length: base.length - length }, (_, i) => `${s.at - i}`));
  return base.filter((l) => !drop.has(l));
}

/**
 * The IMGT framework positions the counting assumes empty (a kappa light chain's FR3 lacks 73, 81 and 82; a heavy
 * chain's FR1 lacks 10 and its FR3 73): every other framework position from the domain's start to 104, and from 118 to
 * its end, holds a residue. Another occupancy (a lambda chain's 81 and 82, a heavy chain's 73, an extra FR2 residue)
 * would put framework residues into a CDR's count, so it is refused, never numbered.
 */
const ImgtFrameworkGaps: Readonly<Record<Chain, readonly number[]>> = { light: [73, 81, 82], heavy: [10, 73] };

function checkFrameworks(chain: Chain, imgt: Labels): void {
  const inCdr = (n: number) => [1, 2, 3].some((c) => n >= ImgtCdrs[c as 1][0] && n <= ImgtCdrs[c as 1][1]);
  const occupied = new Set<number>();
  for (const l of imgt) {
    if (l === null) continue;
    const n = Number.parseInt(l, 10);
    if (inCdr(n)) continue;
    if (`${n}` !== l)
      throw new Error(`${chain}: IMGT ${l} is a framework insertion, which the Kabat counting does not take`);
    occupied.add(n);
  }
  const ns = [...occupied];
  const first = Math.min(...ns);
  const last = Math.max(...ns);
  const gaps = new Set(ImgtFrameworkGaps[chain]);
  for (let p = first; p <= last; p++) {
    if (inCdr(p) || (p > 104 && p < 118)) continue;
    const expected = !gaps.has(p);
    if (occupied.has(p) !== expected)
      throw new Error(
        `${chain}: IMGT ${p} is ${occupied.has(p) ? 'occupied' : 'empty'}; the Kabat counting assumes a ${chain} chain's frameworks with only ${[...gaps].join(', ')} empty`,
      );
  }
}

/** Kabat's labels of a chain from its IMGT labels (`null` where IMGT has none); a framework it does not assume is refused. */
export function kabatFromImgt(chain: Chain, imgt: Labels): (string | null)[] {
  checkFrameworks(chain, imgt);
  const { anchors, stretches } = Kabat[chain];
  const at = anchors.map(([p]) => imgt.indexOf(`${p}`));
  if (at.some((p) => p < 0)) throw new Error(`${chain}: an anchor of IMGT ${anchors.map((a) => a[0])} is missing`);
  const numbered = imgt.flatMap((l, i) => (l === null ? [] : [i]));
  const first = numbered[0];
  const last = numbered[numbered.length - 1];
  const out: (string | null)[] = imgt.map(() => null);
  // before the first cysteine: counted back from it
  for (let i = first; i < at[0]; i++) out[i] = `${anchors[0][1] - (at[0] - i)}`;
  if (out[first] !== null && Number(out[first]) < 1)
    throw new Error(`${chain}: more residues before the cysteine than Kabat numbers`);
  anchors.forEach(([, kabat], k) => {
    out[at[k]] = `${kabat}`;
  });
  stretches.forEach((s, k) => {
    const n = at[k + 1] - at[k] - 1;
    const labels = [...s.before, ...kabatRun(s, n - s.before.length - s.after.length), ...s.after];
    if (labels.length !== n)
      throw new Error(`${chain}: ${n} residues between Kabat ${anchors[k][1]} and ${anchors[k + 1][1]}`);
    labels.forEach((l, j) => {
      out[at[k] + 1 + j] = l;
    });
  });
  // after the last anchor: on to the domain's end
  for (let i = at[3] + 1; i <= last; i++) out[i] = `${anchors[3][1] + (i - at[3])}`;
  return out;
}

/** A CDR by Kabat's definition: its first and last label, by chain (CDR-H1 31–35B). */
export const KabatCdrs: Readonly<Record<Chain, Readonly<Record<1 | 2 | 3, readonly [number, number]>>>> = {
  light: { 1: [24, 34], 2: [50, 56], 3: [89, 97] },
  heavy: { 1: [31, 35], 2: [50, 65], 3: [95, 102] },
};

/** The positions (0-based, half-open) of a CDR in a chain's labels: from its first label to its last one's insertions. */
export function cdrSpan(labels: Labels, first: number, last: number): [start: number, end: number] {
  const num = (l: string | null) => (l === null ? Number.NaN : Number.parseInt(l, 10));
  const start = labels.findIndex((l) => num(l) >= first && num(l) <= last);
  let end = start;
  while (end < labels.length && num(labels[end]) >= first && num(labels[end]) <= last) end++;
  if (start < 0) throw new Error(`no label from ${first} to ${last}`);
  return [start, end];
}
