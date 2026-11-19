/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import { FastaRuntime } from '@molstar/seqx-core/data';
import type { SeqSpecBuilder } from '@molstar/seqx-core/lang';
import type { Derived } from './derived.js';
import { DerivedData } from './derived-data.js';
import { type Design, design } from './design.js';
import { cdrSpan, type Chain, ImgtCdrs, KabatCdrs, kabatFromImgt, type Labels } from './numbering.js';

// The Antibody dataset: trastuzumab's Fab, PDB 1N8Z's light chain (A) and heavy chain (B) read by location from a
// FASTA file, SAbDab's IMGT numbering of their variable domains and Kabat's derived from it each read from a `seq`
// numbering file (a range a residue, its label); the CDRs by both definitions and the design are code. The letters
// and IMGT labels are also held here (`derived-data.ts`, compared with the files by the fixtures test), since the
// CDRs' ranges and the design are written into the specs.

export const Files = {
  fasta: 'fixtures/antibody/1N8Z.fasta',
  imgt: 'fixtures/antibody/1N8Z-imgt.seq',
  kabat: 'fixtures/antibody/1N8Z-kabat.seq',
} as const;

/** The FASTA records' names: the systems of the numbering files. */
export const SystemNames: Readonly<Record<Chain, string>> = { light: '1N8Z_A', heavy: '1N8Z_B' };

/** The chains' lengths as deposited. */
export const Lengths: Readonly<Record<Chain, number>> = { light: 214, heavy: 220 };

export const PdbAttribution = 'wwPDB, PDB 1N8Z (Cho et al. 2003), CC0 1.0';
export const SabdabAttribution =
  'SAbDab (Dunbar et al. 2014; Schneider et al. 2022), IMGT numbering of 1N8Z, CC BY 4.0; Kabat’s derived from it';

/** The data the views hold (the fixture by default). */
export const Data: Derived = DerivedData;

/** The wild type's numberings: IMGT's as published, Kabat's derived. */
export const wildType = (d: Derived = Data) =>
  ({
    light: { letters: d.light.letters, imgt: d.light.imgt, kabat: kabatFromImgt('light', d.light.imgt) },
    heavy: { letters: d.heavy.letters, imgt: d.heavy.imgt, kabat: kabatFromImgt('heavy', d.heavy.imgt) },
  }) as const;

/** The design of the wild type (its edit list the generator's default). */
export const theDesign = (d: Derived = Data): Design => design(d);

/** A chain's letters and numberings, of either instance. */
export interface Numbered {
  readonly letters: string;
  readonly imgt: Labels;
  readonly kabat: Labels;
}

/** What the views read for FASTA (the wild type's letters): its runtime, which the default registry does not hold. */
export const registry = (): SeqRegistryEntry[] => [{ kind: 'format-runtime', runtime: FastaRuntime }];

const ChainLetter: Readonly<Record<Chain, string>> = { light: 'L', heavy: 'H' };

/** A chain's CDRs by a scheme's definition, as ranges rows (`H3`, `CDR-H3 (Kabat)`). */
export function cdrRows(chain: Chain, scheme: 'kabat' | 'imgt', labels: Labels) {
  const spans = ([1, 2, 3] as const).map((c) => {
    const [first, last] = scheme === 'kabat' ? KabatCdrs[chain][c] : ImgtCdrs[c];
    return cdrSpan(labels, first, last);
  });
  const name = (c: number) => `${ChainLetter[chain]}${c}`;
  return {
    start: spans.map((s) => s[0]),
    end: spans.map((s) => s[1]),
    id: spans.map((_, k) => name(k + 1)),
    label: spans.map((_, k) => `CDR-${name(k + 1)} (${scheme === 'kabat' ? 'Kabat' : 'IMGT'})`),
    category: spans.map((_, k) => `CDR${k + 1}`),
  };
}

/**
 * A CDR's address in a scheme's labels, from its first label to its last as the chain has them (a CDR-H1 with
 * insertions ends at `35B`): `/data/heavy/(range :from 95 :to 102 :numbering kabat)`.
 */
export function cdrAddress(chain: Chain, cdr: 1 | 2 | 3, scheme: 'kabat' | 'imgt', labels: Labels): string {
  const rows = cdrRows(chain, scheme, labels);
  const label = (i: number) => {
    const l = labels[i]!;
    // a `.n` label is written quoted (bare, `111.1` is a number)
    return l.includes('.') ? `"${l}"` : l;
  };
  return `/data/${chain}/(range :from ${label(rows.start[cdr - 1])} :to ${label(rows.end[cdr - 1] - 1)} :numbering ${scheme})`;
}

/**
 * The two chains as sequences, each with its Kabat and IMGT numbering: the wild type read from the files (letters from
 * the FASTA, labels from the numbering files), or an instance inline (the design).
 */
export function fab(b: SeqSpecBuilder, instance?: Readonly<Record<Chain, Numbered>>) {
  if (instance) {
    const inline = (chain: Chain) =>
      b.data.sequence(
        chain,
        { alphabet: 'protein', letters: instance[chain].letters, name: `design ${chain}` },
        (s) => {
          s.numbering('kabat', { method: 'labels', labels: [...instance[chain].kabat], scheme: 'kabat' });
          s.numbering('imgt', { method: 'labels', labels: [...instance[chain].imgt], scheme: 'imgt' });
        },
      );
    return { light: inline('light'), heavy: inline('heavy') };
  }
  const fasta = b.data.source('fab', {
    format: 'fasta',
    uri: Files.fasta,
    attribution: PdbAttribution,
    licence: 'CC0-1.0',
  });
  const numbering = (id: 'imgt' | 'kabat') =>
    b.data.source(`${id}-labels`, {
      format: 'seq',
      uri: Files[id],
      attribution: SabdabAttribution,
      licence: 'CC-BY-4.0',
    });
  const imgt = numbering('imgt');
  const kabat = numbering('kabat');
  const read = (chain: Chain) =>
    b.data.sequence(chain, { alphabet: 'protein', source: fasta, name: SystemNames[chain] }, (s) => {
      s.numbering('kabat', { method: 'labels', source: kabat, scheme: 'kabat' });
      s.numbering('imgt', { method: 'labels', source: imgt, scheme: 'imgt' });
    });
  return { light: read('light'), heavy: read('heavy') };
}
