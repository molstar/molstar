/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// AlphaFold DB's PAE JSON as a matrix (`pae.ts` has the format and the provider over it). Nothing imported: the
// model-archive extension's PAE card reads the JSON with it (`@molstar/prism-viewer/sequence/molstar/pae-matrix`)
// without pulling Seq* into the app's start code (sequence.md §2, a runtime chunk).

/**
 * A predicted aligned error matrix: `size` tokens, by scored token (row j: token j's error aligned on each token i, at
 * `j · size + i`; AlphaFold DB's JSON transposed), Å, as float32: a grid's a-major cells, side a scored
 * (`seqspec:model#annotations`).
 */
export interface PaeMatrix {
  readonly kind: 'afdb_pae';
  readonly size: number;
  readonly values: Float32Array;
}
export const PaeMatrix = {
  is: (x: unknown): x is PaeMatrix => (x as PaeMatrix | null)?.kind === 'afdb_pae',
};

/** The matrix of AlphaFold DB's JSON (a list of one object, or the object); a message when it is not one. */
export function readPae(json: unknown): PaeMatrix | string {
  const o = (Array.isArray(json) ? json[0] : json) as Record<string, unknown> | null;
  if (!o || typeof o !== 'object') return 'not an object or a list of one';
  const rows = o.predicted_aligned_error;
  if (Array.isArray(rows)) {
    const n = rows.length;
    const values = new Float32Array(n * n);
    for (let i = 0; i < n; i++) {
      const row = rows[i];
      if (!Array.isArray(row) || row.length !== n) return `row ${i} is not a list of ${n} numbers`;
      for (let j = 0; j < n; j++) values[j * n + i] = typeof row[j] === 'number' ? row[j] : Number.NaN;
    }
    return { kind: 'afdb_pae', size: n, values };
  }
  // The older form: 1-based residue pairs and their distance, every pair listed.
  const [r1, r2, d] = [o.residue1, o.residue2, o.distance];
  if (!Array.isArray(r1) || !Array.isArray(r2) || !Array.isArray(d) || r1.length !== d.length || r2.length !== d.length)
    return 'no predicted_aligned_error, and no residue1, residue2 and distance lists of one length';
  const n = Math.round(Math.sqrt(d.length));
  if (n * n !== d.length) return `${d.length} distances are not N × N`;
  const values = new Float32Array(n * n).fill(Number.NaN);
  for (let k = 0; k < d.length; k++) {
    const [i, j] = [r1[k], r2[k]];
    if (!Number.isInteger(i) || !Number.isInteger(j) || i < 1 || j < 1 || i > n || j > n)
      return `pair ${k} (residue1 ${i}, residue2 ${j}) is not two residues 1..${n}`;
    values[(j - 1) * n + (i - 1)] = typeof d[k] === 'number' ? d[k] : Number.NaN;
  }
  return { kind: 'afdb_pae', size: n, values };
}
