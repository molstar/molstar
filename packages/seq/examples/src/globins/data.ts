/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import { cellsSource, delayedSource, type FormatRuntime } from '@molstar/seqx-core/data';

// The Globins dataset, generated (a stand-in until the Pfam PF00042 seed and its tree are recorded): a family evolved
// from human β-globin (HBB_HUMAN, UniProt P68871, the ancestor) along a random binary tree on core's `PCG`, seeded, so
// every call gives the same rows; its tree as Newick with the branch lengths it evolved along; and a scale variant of
// any size whose rows are made on demand by index from the family's per-column profile, lowercase inserts with an
// `RF`-like column every 25 (EC-38, EC-41), never held all at once by the source.

/** Human β-globin (UniProt P68871), the family's ancestor. */
export const HbbHuman =
  'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH';

const Residues = 'ACDEFGHIKLMNPQRSTVWY';

/** A generated family: its rows (labels, aligned letters, first residue numbers) and its tree. */
export interface Family {
  readonly label: readonly string[];
  readonly letters: readonly string[];
  readonly newick: string;
}

/**
 * The family of `n` rows (48 by default): a random binary tree over `n` leaves, each branch a number of substitutions
 * and, now and then, a run of gaps, from the ancestor down; leaves named `GLB<k>_GEN/<start>-<end>` by their residues
 * (EC-39), the ancestor's row first and named `HBB_HUMAN/1-147`.
 */
export function globinFamily(n = 48, seed = 147): Family {
  const rng = new PCG(seed);
  const rows: string[] = [];
  const parts: string[] = [];
  const evolve = (s: string, steps: number) => {
    const a = s.split('');
    for (let k = 0; k < steps; k++) {
      const j = rng.int() % a.length;
      if (rng.float() < 0.08) {
        const len = 1 + (rng.int() % 4);
        for (let q = j; q < Math.min(a.length, j + len); q++) a[q] = '-';
      } else if (a[j] !== '-') a[j] = Residues[rng.int() % Residues.length];
    }
    return a.join('');
  };
  // builds the subtree over `count` leaves from `seq`: its Newick text, its leaves' rows pushed in order
  const grow = (seq: string, count: number, first: boolean): string => {
    if (count === 1) {
      const k = rows.length;
      rows.push(seq);
      parts.push(first && k === 0 ? 'HBB_HUMAN' : `GLB${k}_GEN`);
      return `#${k}`;
    }
    const left = 1 + (rng.int() % (count - 1));
    const out: string[] = [];
    for (const [c, f] of [
      [left, first],
      [count - left, false],
    ] as const) {
      const steps = f ? 0 : 2 + (rng.int() % 10);
      out.push(`${grow(evolve(seq, steps), c, f)}:${steps}`);
    }
    return `(${out.join(',')})`;
  };
  let newick = `${grow(HbbHuman, n, true)};`;
  const label = rows.map((r, k) => {
    const residues = r.replace(/-/g, '').length;
    return `${parts[k]}/1-${residues}`;
  });
  newick = newick.replace(/#(\d+)/g, (_, k) => label[Number(k)]);
  return { label, letters: rows, newick };
}

/** Each column's residue frequencies and gap rate over a family's rows: what the scale variant samples. */
function profileOf(f: Family): { cdf: Float64Array; gap: Float64Array; columns: number } {
  const C = f.letters[0].length;
  const K = Residues.length;
  const cdf = new Float64Array(C * K);
  const gap = new Float64Array(C);
  for (let j = 0; j < C; j++) {
    let n = 0;
    for (const r of f.letters) {
      const c = r[j];
      if (c === '-') gap[j]++;
      else {
        cdf[j * K + Residues.indexOf(c)]++;
        n++;
      }
    }
    gap[j] /= f.letters.length;
    let acc = 0;
    for (let k = 0; k < K; k++) {
      acc += n ? cdf[j * K + k] / n : 1 / K;
      cdf[j * K + k] = acc;
    }
  }
  return { cdf, gap, columns: C };
}

/**
 * The scale variant's row `i` of `columns` (300 by default): the family's columns first, then columns sampled from its
 * profile; every 25th column an insert (lowercase, gaps mostly); a function of the seed and the index alone.
 */
export function generatedRow(i: number, columns: number, seed: number, profile = defaultProfile(), out?: Uint8Array) {
  const K = Residues.length;
  const row = out ?? new Uint8Array(columns);
  const rng = new PCG((seed * 0x9e3779b1 + i * 2654435761) >>> 0);
  for (let j = 0; j < columns; j++) {
    const pj = j % profile.columns;
    const insert = j % 25 === 24;
    if (rng.float() < (insert ? 0.8 : profile.gap[pj])) {
      row[j] = insert ? 46 : 45; // `.` in an insert column, `-` elsewhere
      continue;
    }
    const u = rng.float();
    let k = 0;
    while (k < K - 1 && profile.cdf[pj * K + k] < u) k++;
    const c = Residues.charCodeAt(k);
    row[j] = insert ? c + 32 : c;
  }
  return row;
}

let profileKept: ReturnType<typeof profileOf> | undefined;
function defaultProfile() {
  return (profileKept ??= profileOf(globinFamily()));
}

/** A row's label in the scale variant: `GEN_<n>/<start>-<end>` (EC-39), its start a function of its index. */
export function generatedLabel(i: number): string {
  const start = 1 + (i % 7);
  return `GEN_${i + 1}/${start}-${start + 140}`;
}

/** A balanced tree over `n` rows of the scale variant (for collapsed clades at scale), leaves by label in data order. */
export function balancedNewick(n: number): string {
  const parts: string[] = [];
  const grow = (a: number, b: number, depth: number): void => {
    if (b - a === 1) {
      parts.push(`${generatedLabel(a)}:1`);
      return;
    }
    const m = (a + b) >> 1;
    parts.push('(');
    grow(a, m, depth + 1);
    parts.push(',');
    grow(m, b, depth + 1);
    parts.push(`)c${a}_${b}:1`);
  };
  grow(0, n, 0);
  return `${parts.join('')};`;
}

/** The scale variant's format: `value: { rows, columns, seed }`, rows made on demand by index (EC-41). */
export const GeneratedFormat = 'globins-generated';

/** Its runtime; each answer after `delayMs` (a story's pacing; 0: none). */
export function generatedRuntime(delayMs = 0): FormatRuntime {
  return {
    name: GeneratedFormat,
    formats: [GeneratedFormat],
    async open(source) {
      const v = (source.value ?? {}) as { rows?: number; columns?: number; seed?: number };
      const rows = v.rows ?? 10_000;
      const columns = v.columns ?? 300;
      const seed = v.seed ?? 1;
      const inner = cellsSource(rows, columns, (start, end) => {
        const codes = new Uint8Array((end - start) * columns);
        const label: string[] = [];
        const first = new Int32Array(end - start);
        for (let i = start; i < end; i++) {
          generatedRow(i, columns, seed, undefined, codes.subarray((i - start) * columns, (i - start + 1) * columns));
          label.push(generatedLabel(i));
          first[i - start] = 1 + (i % 7);
        }
        return { codes, label, first };
      });
      return delayMs > 0 ? delayedSource(inner, { delayMs }) : inner;
    },
  };
}

/** What a view of the scale variant needs: its runtime. */
export function generatedRegistry(options: { readonly delayMs?: number } = {}): readonly SeqRegistryEntry[] {
  return [{ kind: 'format-runtime', runtime: generatedRuntime(options.delayMs ?? 0) }];
}
