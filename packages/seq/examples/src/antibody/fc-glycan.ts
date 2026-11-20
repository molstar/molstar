/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Patch } from '@molstar/core/repr/patch';
import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, identifiers } from '@molstar/seqx-core/lang';

// The Antibody's Fc glycan (examples.md, Antibody; molecules.md): PDB 3AVE, a fucosylated human IgG1 Fc (IGHG1_HUMAN
// P01857 108–330, author numbers 225–447), its two chains each with the N-glycan on Asn297: the branched molecule of
// eight residues (`pdbx_entity_branch_list` in `num` order: NAG NAG BMA MAN NAG MAN NAG FUC, the biantennary G0F) with
// its seven glycosidic links (`pdbx_entity_branch_link`, its own `links`), its WURCS. A complex of the entry's units:
// the Fc chains A and B (a copy group, merged ×2 by default; observed from 234 and 237 to 444), their glycans C and D
// (each its author numbers 1–8) and the zinc E; the complex's links: each glycan's NAG C1 on its chain's Asn297 ND2,
// each chain's two disulfides (C261–C321, C367–C425), the zinc on His268 and His285 of A. The glycan a chip (`whole`)
// or its residues (`residues`, its links drawn as arcs): one patch each way, the editing view's. The Fc's letters and
// chain A's secondary structure from its `seq` file; the rest the entry's, held here and checked against it by a test.

export const title = 'IgG1 Fc and its N297 glycan (3AVE)';

/** The Fc entity's `seq` file (`scripts/record.ts antibody-structure`) and its system there. */
export const FcFile = 'fixtures/antibody/3AVE.seq';
export const FcSystem = '3AVE-1';
export const FcLength = 223;
/** An Fc position's author number: P01857's, less 107 (the entity starts at its 108, author number 225). */
export const fcAuth = (position: number) => position + 225;
/** The glycan's residues in `num` order (1–8, its author numbers too). */
export const GlycanResidues = ['NAG', 'NAG', 'BMA', 'MAN', 'NAG', 'MAN', 'NAG', 'FUC'] as const;
/** Its glycosidic links: the child's C1 (`a`, 0-based) on the parent's oxygen (`b`). */
export const GlycanLinks = [
  { a: 1, b: 0, atom: 'O4' },
  { a: 2, b: 1, atom: 'O4' },
  { a: 3, b: 2, atom: 'O3' },
  { a: 4, b: 3, atom: 'O2' },
  { a: 5, b: 2, atom: 'O6' },
  { a: 6, b: 5, atom: 'O2' },
  { a: 7, b: 0, atom: 'O6' },
] as const;
export const GlycanWurcs =
  'WURCS=2.0/4,8,7/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a1221m-1a_1-5]/1-1-2-3-1-3-1-4/a4-b1_a6-h1_b4-c1_c3-d1_c6-f1_d2-e1_f2-g1';
/** The chains: each one's modeled residues (0-based, end exclusive) and its glycan's unit. */
export const FcChains = [
  { chain: 'A', observed: [9, 220], glycan: 'C' },
  { chain: 'B', observed: [12, 220], glycan: 'D' },
] as const;
/** The glycosylated asparagine (Asn297), its two disulfides, the zinc's histidines: author numbers. */
export const Asn297 = 297;
export const Disulfides = [
  [261, 321],
  [367, 425],
] as const;
export const ZincHistidines = [268, 285] as const;
/** An author number's position on the Fc. */
export const fcPosition = (auth: number) => auth - 225;

export interface FcGlycanOptions {
  /** The glycans drawn as chips (`whole`, absent) or by their residues. */
  readonly glycan?: 'whole' | 'residues';
  readonly copies?: 'merge' | 'each';
}

/** The section's address and the glycan block's path (what the editing view's patch names). */
export const FcSection = '/views/main/fc';
const GlycanBlock = ['views', 'main', 'fc', 'C'];

export function build(options: FcGlycanOptions = {}): Node {
  const b = createSeqSpec({ title });
  const file = b.data.source('3ave-file', {
    format: 'seq',
    uri: FcFile,
    attribution: 'wwPDB, PDB 3AVE (Matsumiya et al. 2007)',
    licence: 'CC0-1.0',
  });
  const fc = b.data.sequence('fc', { alphabet: 'protein', source: file, name: FcSystem });
  const secondary = b.data.ranges('fc-secondary', { target: fc, source: file });
  const glycan = b.data.molecule(
    'glycan',
    {
      class: 'branched',
      name: 'N-glycan (G0F)',
      residues: [...GlycanResidues],
      identifiers: identifiers([['wurcs', GlycanWurcs]]),
    },
    (m) =>
      m.links('links', {
        rows: {
          a_start: GlycanLinks.map((l) => l.a),
          a_end: GlycanLinks.map((l) => l.a + 1),
          b_start: GlycanLinks.map((l) => l.b),
          b_end: GlycanLinks.map((l) => l.b + 1),
          type: GlycanLinks.map(() => 'glycosidic'),
          a_atom: GlycanLinks.map(() => 'C1'),
          b_atom: GlycanLinks.map((l) => l.atom),
        },
      }),
  );
  const zn = b.data.molecule('zn', {
    class: 'ion',
    name: 'zinc',
    residues: ['ZN'],
    atoms: { residue: [0], name: ['ZN'], element: ['ZN'] },
    identifiers: identifiers([['ccd', 'ZN']]),
  });
  const labels = (o: readonly [number, number]) =>
    Array.from({ length: FcLength }, (_, i) => (i < o[0] || i >= o[1] ? null : String(fcAuth(i))));
  b.data.complex('3ave', { name: '3AVE' }, (c) => {
    for (const ch of FcChains)
      c.unit(
        ch.chain,
        { of: fc, chain: ch.chain, auth_chain: ch.chain, observed: { start: [ch.observed[0]], end: [ch.observed[1]] } },
        (u) => u.numbering('auth', { method: 'labels', labels: labels(ch.observed) }),
      );
    for (const ch of FcChains)
      c.unit(ch.glycan, { of: glycan, chain: ch.glycan, auth_chain: ch.chain }, (u) =>
        u.numbering('auth', { method: 'labels', labels: GlycanResidues.map((_, i) => String(i + 1)) }),
      );
    c.unit('E', { of: zn, chain: 'E', auth_chain: 'A' }, (u) =>
      u.numbering('auth', { method: 'labels', labels: ['448'] }),
    );
    const rows = [
      ...FcChains.map((ch) => ({
        a: ch.chain,
        ap: fcPosition(Asn297),
        aa: 'ND2',
        b: ch.glycan,
        bp: 0,
        ba: 'C1',
        type: 'glycosidic',
      })),
      ...FcChains.flatMap((ch) =>
        Disulfides.map(([x, y]) => ({
          a: ch.chain,
          ap: fcPosition(x),
          aa: 'SG',
          b: ch.chain,
          bp: fcPosition(y),
          ba: 'SG',
          type: 'disulfide',
        })),
      ),
      ...ZincHistidines.map((h) => ({ a: 'A', ap: fcPosition(h), aa: 'ND1', b: 'E', bp: 0, ba: 'ZN', type: 'metal' })),
    ];
    c.links('links', {
      rows: {
        a_unit: rows.map((r) => r.a),
        a_start: rows.map((r) => r.ap),
        a_end: rows.map((r) => r.ap + 1),
        b_unit: rows.map((r) => r.b),
        b_start: rows.map((r) => r.bp),
        b_end: rows.map((r) => r.bp + 1),
        type: rows.map((r) => r.type),
        a_atom: rows.map((r) => r.aa),
        b_atom: rows.map((r) => r.ba),
      },
    });
  });
  b.views.view('main', (v) => {
    const s = v.section('fc', { of: '/data/3ave' as never, ...(options.copies ? { copies: options.copies } : {}) });
    // the glycans' block (C, the merged copies' first unit) overridden alone: its form, the rest derived
    s.block('C', { form: options.glycan ?? 'whole' });
    s.track('links', { label: 'Links' }).layer('links', { type: 'pairs', data: '/data/3ave/links' as never });
    s.track('glycan', { label: 'Glycan links' }).layer('glycan', {
      type: 'pairs',
      data: '/data/glycan/links' as never,
    });
    s.track('secondary', { label: 'Secondary structure' }).layer('secondary', { type: 'ranges', data: secondary });
    s.track('letters', { label: 'Fc' }).layer('letters', { type: 'letters' });
  });
  return b.build();
}

/** The glycans by their residues, their glycosidic links drawn. */
export const residues = () => build({ glycan: 'residues' });

/** The editing view's toggle: the glycans as chips or by their residues, one patch. */
export function glycanFormPatch(form: 'whole' | 'residues'): Patch {
  return { ops: [{ op: 'set', at: GlycanBlock, key: 'form', value: form }] };
}
