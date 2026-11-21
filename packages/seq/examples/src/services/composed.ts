/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, type SeqSpecBuilder, type ViewBuilder } from '@molstar/seqx-core/lang';
import { nodeAddress } from '@molstar/seqx-core/resolve';
import { plddtColor } from '../complex/confidence.js';
import { uniprotEntry } from '../files/uniprot.js';
import { AlphaFoldDb, Pdbe, serviceUrl } from './endpoints.js';
import { addSiftsMaps, chainsOf, type SiftsMapping } from './sifts.js';

// An application composing services (examples.md, Showcases): p53 (P04637) from UniProt, its per-residue confidence from
// AlphaFold DB, PDB 2OCJ from PDBe, and PDBe's SIFTS mappings between them, each by the location its service answers at
// (`endpoints.ts`): the application's assets port answers from recorded fixtures or from the services themselves, the
// spec is the same. UniProt's and AlphaFold DB's answers are read by the spec's own sources (records, a transform each);
// SIFTS's is read by the application (`sifts.ts`) into one `structure` map a chain, onto the unit the mmCIF source
// describes (a described end: pending until the structure has arrived), so the view is first drawn without the maps'
// coverage and the structure section, and drawn again once SIFTS has answered. `compose` writes what every scenario of
// the showcase shares (the predictions', the addressing's).

export const title = 'An application composing UniProt, PDBe and AlphaFold DB';

export const Accession = 'P04637';
export const Entry = '2OCJ';

/** The addresses the application and its adapter name. */
export const ServiceAddresses = {
  protein: '/data/protein',
  entry: '/data/entry',
  confidence: '/data/plddt',
  structure: '/data/2ocj',
  complex: nodeAddress('/data/2ocj', 'complex'),
} as const;

/**
 * A chain's unit in the complex the structure source describes (`<source>/(node :key complex)/<chain>`): where SIFTS's
 * maps end (a described end, pending until the description arrives) and what the structure link translates to.
 */
export const unitAddress = (chain: string) => `${ServiceAddresses.complex}/${chain}`;

/** The spec's sequence of a chain of the structure source (named by the source: `source`, the unit's id as `name`). */
export const chainAddress = (chain: string) => `/data/2ocj-${chain}`;

export interface ComposedOptions {
  /** SIFTS's answer for the entry; absent while it is on its way (the view without it). */
  readonly sifts?: readonly SiftsMapping[];
}

/** What `compose` wrote: the addresses its nodes have. */
export interface Composed {
  readonly protein: string;
  readonly plddt: string;
  /** The ranges of the entry's chains on the protein (SIFTS's), absent until it has answered. */
  readonly coverage?: string;
  /** The data of a layer of one UniProt feature type (`Domain`, `Region`). */
  readonly features: (type: string) => string;
  /** Whether SIFTS has answered: its maps and the structure's sequences are written. */
  readonly structure: boolean;
  /** The B-factors of chain A (PDBe's mmCIF), when the structure is written. */
  readonly bfactor?: string;
}

/** The nodes every scenario shares: the three services' data and, when SIFTS has answered, its maps. */
export function compose(b: SeqSpecBuilder, options: ComposedOptions = {}): Composed {
  const uniprot = uniprotEntry(b, serviceUrl.uniprot(Accession), { accession: Accession });
  const confidence = b.data.source('plddt', {
    format: 'json',
    uri: serviceUrl.confidence(Accession),
    fields: { position: 'residueNumber', value: 'confidenceScore' },
    system: Accession,
    attribution: `${AlphaFoldDb.attribution} AF-${Accession}-F1 v6`,
    licence: AlphaFoldDb.licence,
  });
  const plddt = b.data.values('confidence', { target: uniprot.protein, source: confidence });
  const pdb = b.data.source('2ocj', {
    format: 'mmcif',
    uri: serviceUrl.entry(Entry),
    attribution: `${Pdbe.attribution}, PDB ${Entry}`,
    licence: 'CC0-1.0',
  });
  const sifts = options.sifts ?? [];
  const rows = sifts.filter((m) => m.accession === Accession);
  const coverage = rows.length
    ? b.data.ranges('pdb-coverage', {
        target: uniprot.protein,
        rows: {
          start: rows.map((m) => m.unpStart - 1),
          end: rows.map((m) => m.unpEnd),
          id: rows.map((m) => `${m.entry} ${m.chain}`),
          label: rows.map((m) => `${m.entry} chain ${m.chain} (${m.unpStart}–${m.unpEnd})`),
        },
      })
    : undefined;
  for (const chain of new Set(rows.map((m) => m.chain)))
    b.data.sequence(`2ocj-${chain}`, { alphabet: 'protein', source: pdb, name: chain });
  addSiftsMaps(b, Entry, rows, { [Accession]: uniprot.protein }, unitAddress);
  const bfactor = rows.length
    ? b.data.values('bfactor-a', { target: chainAddress('A') as never, source: pdb, field: 'b_factor' })
    : undefined;
  return {
    protein: uniprot.protein,
    plddt,
    ...(coverage ? { coverage } : {}),
    features: uniprot.data,
    structure: rows.length > 0,
    ...(bfactor ? { bfactor } : {}),
  };
}

/**
 * The section of chain A of the structure, once SIFTS has answered: PDBe's mmCIF serves a unit's B-factors (a residue
 * without atoms has none: the gaps at the ends are the crystal's), never its letters.
 */
export function structureSection(v: ViewBuilder, c: Composed) {
  if (!c.structure || !c.bfactor) return;
  v.section('structure', { of: chainAddress('A') as never }, (s) => {
    s.track('ruler', { label: `${Entry} chain A (PDBe)` }).layer('ruler', { type: 'ruler' });
    s.track('bfactor', { label: 'B-factor (Å²)' }).layer('bfactor', {
      type: 'values',
      data: c.bfactor!,
      domain: [0, 80],
    });
  });
}

/** The spec: three services' data, SIFTS's maps when its answer is here. */
export function build(options: ComposedOptions = {}): Node {
  const b = createSeqSpec({ title });
  const c = compose(b, options);
  b.views.view('main', (v) => {
    v.section('protein', { of: c.protein }, (s) => {
      s.track('ruler', { label: 'P04637' }).layer('ruler', { type: 'ruler' });
      s.track('letters', { label: 'Sequence (UniProt)' }).layer('letters', { type: 'letters' });
      s.track('dna', { label: 'DNA binding (UniProt)' }).layer('dna', {
        type: 'ranges',
        data: c.features('DNA binding'),
      });
      s.track('plddt', { label: 'pLDDT (AlphaFold DB)' }).layer('plddt', {
        type: 'values',
        data: c.plddt,
        color: plddtColor(),
      });
      if (c.coverage)
        s.track('coverage', { label: `${Entry} chains (SIFTS)` }).layer('coverage', {
          type: 'ranges',
          data: c.coverage,
        });
    });
    structureSection(v, c);
  });
  return b.build();
}

/** The chains SIFTS maps onto the accession: `A, B, C, D`. */
export const chainsLine = (sifts: readonly SiftsMapping[]) => chainsOf(sifts, Accession).join(', ');
