/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, type SeqSpecBuilder, type ViewBuilder } from '@molstar/seqx-core/lang';
import { nodeAddress } from '@molstar/seqx-core/resolve';
import { plddtColor } from '../complex/confidence.js';
import { afdbConfidence } from '../complex/data.js';
import { uniprotEntry } from '../files/uniprot.js';
import { Pdbe, serviceUrl } from './endpoints.js';
import { addSiftsMaps, chainsOf, isOneToOne, type SiftsMapping, siftsSegments } from './sifts.js';

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
  /** The unit the structure section shows (`struct_asym_id`), and the static map the section follows through. */
  readonly asym?: string;
  readonly via?: string;
}

/** The nodes every scenario shares: the three services' data and, when SIFTS has answered, its maps. */
export function compose(b: SeqSpecBuilder, options: ComposedOptions = {}): Composed {
  const uniprot = uniprotEntry(b, serviceUrl.uniprot(Accession), { accession: Accession });
  const plddt = afdbConfidence(b, {
    source: 'plddt',
    values: 'confidence',
    accession: Accession,
    uri: serviceUrl.confidence(Accession),
    target: uniprot.protein,
  });
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
  // the units are `struct_asym_id`'s (what the mmCIF source calls them), never the author's `chain_id`
  const asyms = [...new Set(siftsSegments(rows).map((s) => s.asym))];
  // each unit is numbered by the author's numbers of its mapped run (SIFTS's, else UniProt's: 2OCJ's author numbers are
  // UniProt's), so its ruler counts 94–312 under P04637's, never the entity's own 1–219
  for (const asym of asyms) {
    const runs = rows.filter((m) => m.asym === asym && isOneToOne(m));
    b.data.sequence(`2ocj-${asym}`, { alphabet: 'protein', source: pdb, name: asym }, (s) => {
      s.numbering('auth', {
        method: 'segments',
        scheme: 'auth',
        segments: {
          start: runs.map((m) => m.start - 1),
          end: runs.map((m) => m.end),
          first: runs.map((m) => m.authStart ?? m.unpStart),
        },
      } as never);
    });
  }
  // two sets of maps, one a chain: onto the units the structure link translates to (a described end, pending until the
  // structure has arrived), and onto the spec's own sequences of the units, which the sections show and the hover and the
  // selection coming back from the structure are projected through (the static model holds those, not the units)
  addSiftsMaps(b, Entry, rows, { [Accession]: uniprot.protein }, unitAddress);
  const onSequences = addSiftsMaps(b, Entry, rows, { [Accession]: uniprot.protein }, chainAddress, '-seq');
  const asym = asyms[0];
  const bfactor = asym
    ? b.data.values('bfactor', { target: chainAddress(asym) as never, source: pdb, field: 'b_factor' })
    : undefined;
  return {
    protein: uniprot.protein,
    plddt,
    ...(coverage ? { coverage } : {}),
    features: uniprot.data,
    structure: rows.some(isOneToOne),
    ...(bfactor && asym ? { bfactor, asym, via: [...onSequences.values()][0] } : {}),
  };
}

/**
 * The section of chain A of the structure, once SIFTS has answered: PDBe's mmCIF serves a unit's B-factors (a residue
 * without atoms has none: the gaps at the ends are the crystal's), never its letters.
 */
export function structureSection(v: ViewBuilder, c: Composed, follow: string) {
  if (!c.structure || !c.bfactor || !c.asym) return;
  // the section follows `follow` (the protein's, or the lab's) through SIFTS's map: the B-factors read under the residues
  // of P04637 they are paired with, the ruler counting the author's numbers, which are P04637's
  v.section(
    'structure',
    {
      of: chainAddress(c.asym) as never,
      follow: follow as never,
      ...(c.via ? { via: c.via as never } : {}),
    },
    (s) => {
      s.track('ruler', { label: `${Entry} chain ${c.asym} (PDBe, author numbers)` }).layer('ruler', {
        type: 'ruler',
        numbering: 'auth',
      });
      s.track('bfactor', { label: 'B-factor (Å²)' }).layer('bfactor', {
        type: 'values',
        data: c.bfactor!,
        domain: [0, 80],
      });
    },
  );
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
    structureSection(v, c, '/views/main/protein');
  });
  return b.build();
}

/** The chains SIFTS maps onto the accession: `A, B, C, D`. */
export const chainsLine = (sifts: readonly SiftsMapping[]) => chainsOf(sifts, Accession).join(', ');
