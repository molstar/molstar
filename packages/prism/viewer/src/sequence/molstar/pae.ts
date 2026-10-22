/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// AlphaFold DB's predicted aligned error (seqspec:sources#formats-and-services, seqspec:libraries#the-mol-adapter): the
// `afdb_pae` format reads AlphaFold DB's PAE JSON (`predicted_aligned_error`, an N × N list; its
// `max_predicted_aligned_error` is not read; the older `residue1`, `residue2`, `distance` lists alike) into an N × N
// float32 matrix stored by scored token (transposed once as read), and `molstar:pae` reads it into a `grid` over its
// complex's polymer units in the complex's order (ligands are no tokens of it), a row per pair of units, each a
// Float32Array the grid keeps as given (a single chain's the matrix itself, no copy): side a the scored residue, side
// b the residue aligned on, so a dot plot draws it as AlphaFold DB does (scored across, aligned down). No Mol* here:
// the matrix's tokens are the units' residues in order.

import {
  type FormatEntry,
  issueOf,
  type ProviderEntry,
  RequiredField,
  ref,
  SimpleParamsSchema,
  type ThemeRef,
} from '@molstar/seq-core/format';
import { KeyKind } from '@molstar/seq-core/model';
import { PaeMatrix, readPae } from './pae-matrix.js';

export { PaeMatrix, readPae };

/** AlphaFold DB's PAE plot colors: dark green at 0 Å to white at 30 Å and beyond. */
export const PaeScale: ThemeRef = {
  theme: 'field_palette',
  params: {
    field: 'value',
    palette: {
      kind: 'continuous',
      mode: 'absolute',
      colors: [
        ['#0e5326', 0],
        ['#3b9a54', 10],
        ['#a5d8a7', 20],
        ['#f4fbf2', 30],
      ],
      overflow_color: '#f4fbf2',
    },
  },
};

export const PaeFormat: FormatEntry = {
  name: 'afdb_pae',
  extensions: [],
  gives: ['object'],
  reads: [],
  read(bytes, _request, ctx) {
    let json: unknown;
    try {
      json = JSON.parse(new TextDecoder().decode(bytes));
    } catch {
      json = null;
    }
    const m = readPae(json);
    if (typeof m !== 'string') return { kind: 'object', value: m };
    ctx.report(issueOf(ctx.registries, 'svs.source.format', '', `AlphaFold DB's PAE: ${m}`));
    return { kind: 'object', value: null };
  },
};

export const PaeProvider: ProviderEntry = {
  name: 'molstar:pae',
  options: SimpleParamsSchema({
    of: RequiredField(ref('complex'), 'the complex of a predicted model: its polymer units, in order, are the tokens'),
    source: RequiredField(ref('source'), "AlphaFold DB's PAE (service afdb, resource pae)"),
  }),
  geometry: 'grid',
  defaultColor: PaeScale,
  async compute(ctx, options) {
    const o = options as { of: string; source: string };
    const parsed = ctx.source(o.source);
    const m = parsed?.kind === 'object' && PaeMatrix.is(parsed.value) ? parsed.value : null;
    const set = ctx.set;
    const t = { kind: 0, id: 0, row: -1 };
    if (!m || !set) return null;
    if (!set.byName(o.of, t) || t.kind !== KeyKind.complex) {
      const message = `"${o.of}" is not a complex: give the predicted model's complex, whose polymer units are the PAE's tokens`;
      ctx.report(issueOf(ctx.registries, 'svs.provider.unavailable', '', message));
      return null;
    }
    // The complex's polymer units, each over its range, in the complex's order.
    const units = Array.from(set.complexes[t.id]?.units ?? []).filter(
      (u) => !set.molecules[set.units.molecule[u]]?.element,
    );
    const start = units.map((u) => set.units.rangeStart[u]);
    const length = units.map((u, k) => set.units.rangeStop[u] - start[k]);
    const total = length.reduce((a, b) => a + b, 0);
    if (total !== m.size) {
      const message = `the PAE has ${m.size} tokens, the complex's polymer units ${total} residues`;
      ctx.report(issueOf(ctx.registries, 'svs.source.format', '', message));
      return null;
    }
    const cols = { target: [] as string[], start: [] as number[], stop: [] as number[], value: [] as Float32Array[] };
    const b = { b_target: [] as string[], b_start: [] as number[], b_stop: [] as number[] };
    // Each unit's first token: a running sum.
    const offset: number[] = [];
    for (let k = 0, sum = 0; k < length.length; sum += length[k], k++) offset.push(sum);
    for (let p = 0; p < units.length; p++)
      for (let q = 0; q < units.length; q++) {
        // Scored residue i of unit p's row holds its errors aligned on unit q's residues: one run of the matrix's row.
        const one = units.length === 1;
        const block = one ? m.values : new Float32Array(length[p] * length[q]);
        for (let i = 0; !one && i < length[p]; i++) {
          const from = (offset[p] + i) * m.size + offset[q];
          block.set(m.values.subarray(from, from + length[q]), i * length[q]);
        }
        cols.target.push(set.strings.get(set.units.key[units[p]]));
        cols.start.push(start[p]);
        cols.stop.push(start[p] + length[p]);
        b.b_target.push(set.strings.get(set.units.key[units[q]]));
        b.b_start.push(start[q]);
        b.b_stop.push(start[q] + length[q]);
        cols.value.push(block);
      }
    return { key: '', geometry: 'grid', semanticType: 'pae', label: 'PAE', columns: { ...cols, ...b } };
  },
};
