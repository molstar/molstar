/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// I1's shared fixtures: the scenes, small documents with ligands and poses made in the document (no recorded files),
// for the format half and the viewer half alike, each a fingerprint of a different kind of system: a unit, a sequence,
// an alignment's columns, a block reached through a mapping; many poses and kinds; mixed data and the absent property;
// null categories. Frames are built in the tests that import this afresh (heap guards).

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import type { SeqViewerRegistryEntry } from '../../src/draw/api.js';
import {
  type Fingerprint,
  FingerprintProp,
  FingerprintProvider,
  fingerprint,
  interactionColors,
  interactionGlyphs,
} from '../../src/interactions/index.js';

type Builder = ReturnType<typeof createSeqViewSpec>;

export const KIND_NAMES = [
  'hbond',
  'pi_stacking',
  'cation_pi',
  'hydrophobic',
  'salt_bridge',
  'halogen',
  'metal',
  'water_bridge',
] as const;

/** One pose's contacts: a residue index of its protein and a kind (null: a row with no category). */
export type Contacts = readonly (readonly [position: number, kind: string | null])[];

const RESIDUES = 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQYNCDPGWAHM';
export const residues = (n: number, from = 0) =>
  Array.from({ length: n }, (_, i) => RESIDUES[(from + i) % RESIDUES.length]).join('');

/** `n` ligand instances `${prefix}${k}` of one ligand (four heavy atoms). */
function poses(b: Builder, prefix: string, n: number) {
  const ligand = b.data.ligand({ ref: `${prefix}lig`, value: ['C', 'N', 'O', 'C'] });
  return Array.from({ length: n }, (_, k) =>
    b.data.unit({ ref: `${prefix}${k}`, ligand, label: `${prefix} pose ${k + 1}` } as never),
  );
}

/** The pair rows of `contacts[k]` (pose k) against a protein unit `target`. */
function pairs(b: Builder, ref: string, target: string, prefix: string, contacts: readonly Contacts[], more = {}) {
  const rows = contacts.flatMap((cs, k) => cs.map(([position, kind], a) => [k, position, kind, a % 4] as const));
  return b.data.annotation({
    ref,
    geometry: 'pairs',
    value: {
      target: rows.map(() => target),
      start: rows.map((r) => r[1]),
      b_target: rows.map((r) => `${prefix}${r[0]}`),
      b_start: rows.map((r) => r[3]),
      category: rows.map((r) => r[2]),
    },
    ...more,
  } as never);
}

export interface Scene {
  readonly doc: SeqViewSpecDoc;
  /** The section to open. */
  readonly section: string;
  /** The track of the fingerprint layer. */
  readonly track: string;
}

const LAYER = { type: 'interactions:fingerprint', fallback: 'matrix' } as const;

/** Poses against residues 3, 7, 12 and 20: the first two hold hydrogen bonds at 7 and 12, the third lost the one at 7. */
export const BASIC: readonly Contacts[] = [
  [
    [3, 'hydrophobic'],
    [7, 'hbond'],
    [12, 'hbond'],
    [20, 'pi_stacking'],
  ],
  [
    [3, 'hydrophobic'],
    [7, 'hbond'],
    [12, 'hbond'],
    [20, 'pi_stacking'],
    [20, 'hydrophobic'],
  ],
  [
    [3, 'hydrophobic'],
    [12, 'hbond'],
    [20, 'cation_pi'],
  ],
  [
    [7, 'weird'],
    [12, null],
  ],
];

const ref = (n: { ref?: string } | string) => (typeof n === 'string' ? n : String(n.ref));

/** One fingerprint layer's track after a ruler, in a view `v`; `hint`: the track's `height` in rows. */
function view(b: Builder, block: (v: ReturnType<Builder['view']>) => void, data: unknown, hint?: number, map = false) {
  const v = b.view({ ref: 'v', label: 'V' });
  block(v);
  v.track({ ref: 'v.ruler', label: '' }).ruler();
  v.track({ ref: 'v.fp', label: 'FP', ...(hint === undefined ? {} : { height: hint }) } as never).layer({
    ...LAYER,
    data,
  } as never);
  if (map) v.track({ ref: 'v.map', label: 'Map' }).layer({ type: 'seqv:map', fallback: 'omit' } as never);
}

/** The pairs as arcs on the contacts of a view, in the extension's themes (the kernel fallbacks draw them without). */
function arcs(b: Builder, c: unknown) {
  const v = b.view({ ref: 'w', label: 'W' });
  v.track({ ref: 'w.arcs', label: 'Arcs' }).pairs({
    data: c,
    color: interactionColors(),
    glyph: interactionGlyphs(),
  } as never);
  return v;
}

/** The fingerprint over a unit: its columns are the unit's own positions, the block the unit. */
export function unitScene(
  o: { contacts?: readonly Contacts[]; hint?: number; arcs?: boolean; map?: boolean } = {},
): Scene {
  const contacts = o.contacts ?? BASIC;
  const b = createSeqViewSpec({ title: 'Unit system' });
  const seq = b.data.sequence({ ref: 's', label: 'S', value: residues(40) });
  const u = b.data.unit({ ref: 'P', sequence: seq, label: 'P' } as never);
  poses(b, 'L', contacts.length);
  const c = pairs(b, 'c', 'P', 'L', contacts);
  const fp = fingerprint(b.data, { ref: 'fp', pairs: [c], ligands: contacts.map((_, k) => `L${k}`), system: u });
  view(b, (v) => v.block({ unit: u }), fp, o.hint, o.map);
  if (o.arcs) arcs(b, c).block({ unit: u });
  return { doc: b.getDocument(), section: 'v.main', track: 'v.fp' };
}

/** The system a sequence, the block a unit of it; 20 poses with 8 kinds each, every kind at its own residues. */
export function manyScene(poses_ = 20, kinds = 8, hint?: number): Scene {
  const contacts = Array.from({ length: poses_ }, (_, k) =>
    Array.from({ length: kinds }, (_, t): readonly [number, string] => [2 + ((k + 3 * t) % 30), KIND_NAMES[t]]),
  );
  const b = createSeqViewSpec({ title: 'Many poses' });
  const seq = b.data.sequence({ ref: 's', label: 'S', value: residues(40) });
  const u = b.data.unit({ ref: 'P', sequence: seq, label: 'P' } as never);
  poses(b, 'L', poses_);
  const c = pairs(b, 'c', 'P', 'L', contacts);
  const fp = fingerprint(b.data, { ref: 'fp', pairs: [c], ligands: contacts.map((_, k) => `L${k}`), system: seq });
  view(b, (v) => v.block({ unit: u }), fp, hint);
  return { doc: b.getDocument(), section: 'v.main', track: 'v.fp' };
}

/** The system an alignment's columns (hover kind `column`); two rows, the second with a gap at column 5. */
export function alignmentScene(): Scene {
  const b = createSeqViewSpec({ title: 'Alignment system' });
  const r0 = residues(30);
  const r1 = `${r0.slice(0, 5)}-${r0.slice(5)}`;
  const s0 = b.data.sequence({ ref: 's0', label: 'S0', value: r0 });
  const s1 = b.data.sequence({ ref: 's1', label: 'S1', value: r0 });
  const u0 = b.data.unit({ ref: 'P0', sequence: s0, label: 'P0' } as never);
  const u1 = b.data.unit({ ref: 'P1', sequence: s1, label: 'P1' } as never);
  const aln = b.data.alignment({ ref: 'aln', value: [r0.padEnd(31, '-'), r1], labels: ['R0', 'R1'] } as never);
  aln.bind(0, u0, { method: 'align' });
  aln.bind(1, u1, { method: 'align' });
  poses(b, 'L', 3);
  // Pose 0 and 1 contact the first chain, pose 2 the second: the same residues land in different columns.
  const c0 = pairs(b, 'c0', 'P0', 'L', [
    [
      [3, 'hbond'],
      [10, 'hydrophobic'],
    ],
    [[3, 'hbond']],
    [],
  ]);
  const c1 = pairs(b, 'c1', 'P1', 'L', [
    [],
    [],
    [
      [10, 'hbond'],
      [15, 'pi_stacking'],
    ],
  ]);
  const fp = fingerprint(b.data, { ref: 'fp', pairs: [c0, c1], ligands: ['L0', 'L1', 'L2'], system: aln });
  view(b, (v) => v.block({ alignment: aln }), fp);
  return { doc: b.getDocument(), section: 'v.main', track: 'v.fp' };
}

/** The system a sequence; the block another sequence mapped onto it with an offset of 5 (a block reached through a mapping). */
export function mappingScene(): Scene {
  const b = createSeqViewSpec({ title: 'Through a mapping' });
  const s1 = b.data.sequence({ ref: 's1', label: 'S1', value: residues(40) });
  const s2 = b.data.sequence({ ref: 's2', label: 'S2', value: residues(30, 5) });
  const u = b.data.unit({ ref: 'P', sequence: s1, label: 'P' } as never);
  b.data.mapping({ ref: 'm', from: s2, to: s1, method: 'segments', segments: [[0, 5, 30]] } as never);
  poses(b, 'L', 3);
  // s1 positions 8, 12, 20 and 30 are s2's 3, 7, 15 and 25.
  const c = pairs(b, 'c', 'P', 'L', [
    [
      [8, 'hbond'],
      [20, 'hydrophobic'],
    ],
    [
      [8, 'hbond'],
      [12, 'pi_stacking'],
    ],
    [[30, 'hbond']],
  ]);
  const fp = fingerprint(b.data, { ref: 'fp', pairs: [c], ligands: ['L0', 'L1', 'L2'], system: s1 });
  view(b, (v) => v.block({ sequence: s2 }), fp);
  void u;
  return { doc: b.getDocument(), section: 'v.main', track: 'v.fp' };
}

/**
 * A provider of a fingerprint that carries fewer rows than its matrix has states (a host's annotation with the runtime
 * property attached): `ix:short` is the extension's provider with the property cut to `rows` ligands.
 */
export const ShortFingerprint: SeqViewerRegistryEntry = {
  svs: {
    providers: [
      {
        ...FingerprintProvider,
        name: 'ix:short',
        async compute(ctx, options) {
          const a = await FingerprintProvider.compute(ctx, options);
          const fp = a?.props?.[FingerprintProp] as Fingerprint | undefined;
          if (!a || !fp) return a;
          const C = fp.columns.length;
          const cut: Fingerprint = { ...fp, rows: fp.rows.slice(0, 3), cells: fp.cells.slice(0, 3 * C) };
          return { ...a, props: { [FingerprintProp]: cut } };
        },
      },
    ],
  },
};

/**
 * Three blocks on one layer, five poses each (a layer's matrices share their states): the first a fingerprint, the
 * second a fingerprint cut to its first three rows (`ShortFingerprint`), the third a plain matrix.
 */
export function mixedScene(): Scene {
  const b = createSeqViewSpec({ title: 'Mixed data' });
  const block = (key: string, from: number) => {
    const s = b.data.sequence({ ref: `s${key}`, label: key, value: residues(30, from) });
    return b.data.unit({ ref: `P${key}`, sequence: s, label: key } as never);
  };
  const ua = block('A', 0);
  const ub = block('B', 5);
  const uc = block('C', 9);
  poses(b, 'a', 5);
  const states = [0, 1, 2, 3, 4].map((k) => `a${k}`);
  const ca = pairs(b, 'ca', 'PA', 'a', [
    [[3, 'hbond']],
    [
      [3, 'hbond'],
      [9, 'hydrophobic'],
    ],
    [[9, 'hydrophobic']],
    [[3, 'pi_stacking']],
    [[9, 'hbond']],
  ]);
  const cb = pairs(b, 'cb', 'PB', 'a', [
    [
      [4, 'hbond'],
      [9, 'salt_bridge'],
    ],
    [[4, 'hydrophobic']],
    [[9, 'salt_bridge']],
    [[4, 'pi_stacking']],
    [[9, 'hbond']],
  ]);
  const fa = fingerprint(b.data, { ref: 'fa', pairs: [ca], ligands: states, system: ua });
  const fb = b.data.computed({
    ref: 'fb',
    provider: 'ix:short',
    options: { pairs: [cb], ligands: states, system: ub },
  } as never);
  const plain = plainMatrix(b, 'plain', uc, [3, 4, 5, 9], states);
  view(b, (v) => v.block({ unit: ua }).block({ unit: ub }).block({ unit: uc }), [fa, fb, plain]);
  return { doc: b.getDocument(), section: 'v.main', track: 'v.fp' };
}

/** A matrix annotation no fingerprint provider computed: four states, a value (the kind count) at each listed position. */
function plainMatrix(b: Builder, name: string, target: unknown, positions: number[], states = ['w', 'x', 'y', 'z']) {
  const index = positions;
  const value = positions.map((p) => states.map((_, q) => ((p + q) % 3 === 0 ? null : 1 + ((p + q) % 3))));
  return b.data.annotation({ ref: name, geometry: 'matrix', target, states, value: { index, value } } as never);
}

/** The property absent: the layer's whole data is a plain matrix. */
export function absentScene(hint?: number): Scene {
  const b = createSeqViewSpec({ title: 'No fingerprint' });
  const s = b.data.sequence({ ref: 's', label: 'S', value: residues(30) });
  const u = b.data.unit({ ref: 'P', sequence: s, label: 'P' } as never);
  const plain = plainMatrix(b, 'plain', u, [3, 4, 5, 9, 14]);
  view(b, (v) => v.block({ unit: u }), plain, hint);
  return { doc: b.getDocument(), section: 'v.main', track: 'v.fp' };
}

/** A long sequence (20,000 positions) with contacted residues spread over it: many cells under a device pixel. */
export function wideScene(length = 20_000, system: 'unit' | 'sequence' = 'unit'): Scene {
  const contacts = Array.from({ length: 6 }, (_, k) =>
    Array.from({ length: 150 }, (_, i): readonly [number, string] => [
      Math.floor(((i * 131 + k * 17) * length) / 20_000) % length,
      KIND_NAMES[(i + k) % 5],
    ]),
  );
  const b = createSeqViewSpec({ title: 'Wide' });
  const seq = b.data.sequence({ ref: 's', label: 'S', value: residues(length) });
  const u = b.data.unit({ ref: 'P', sequence: seq, label: 'P' } as never);
  poses(b, 'L', contacts.length);
  const c = pairs(b, 'c', 'P', 'L', contacts);
  const ligands = contacts.map((_, k) => `L${k}`);
  const fp = fingerprint(b.data, { ref: 'fp', pairs: [c], ligands, system: system === 'unit' ? u : seq });
  view(b, (v) => v.block({ unit: u }), fp);
  return { doc: b.getDocument(), section: 'v.main', track: 'v.fp' };
}

void ref;
