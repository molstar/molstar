/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { categoryColumn } from '@molstar/core/repr/tree';
import { encodeSeq, inlineSeq, type SeqData } from '@molstar/seqx-core/data';
import { cifRows } from '../files/cif.js';
import type { Chain, Labels } from './numbering.js';

// What the record script keeps (`scripts/record.ts antibody`): PDB 1N8Z's light and heavy chains (A, B) as FASTA, their
// letters as the entry's `_entity_poly` has them; SAbDab's IMGT numbering of each variable domain (its numbered
// segment: a label a residue, `-` an IMGT position the domain has no residue at) as a label a position of the chain;
// and a numbering's `seq` file, a range a labeled position with its label. Here, in the sources, so a test holds the
// rules; the script imports them.

export const Entry = '1N8Z';

/** The chains taken (HER2, chain C, is not), their FASTA names (the systems of the numbering files). */
export const FabChains: Readonly<Record<Chain, { readonly asym: string; readonly name: string }>> = {
  light: { asym: 'A', name: '1N8Z_A' },
  heavy: { asym: 'B', name: '1N8Z_B' },
};

/** The two chains' letters and the entities' descriptions, as the entry's `_entity_poly` and `_entity` give them. */
export function fabChains(cif: string): Record<Chain, { readonly letters: string; readonly description: string }> {
  const entities = new Map(cifRows(cif, '_entity').map((r) => [r.id, r.pdbx_description]));
  const poly = cifRows(cif, '_entity_poly');
  const of = (asym: string) => {
    const r = poly.find((p) => p.pdbx_strand_id.split(',').includes(asym));
    if (!r) throw new Error(`${Entry}: no polymer of chain ${asym}`);
    return {
      letters: r.pdbx_seq_one_letter_code_can.replaceAll(/\s/g, ''),
      description: entities.get(r.entity_id) ?? '',
    };
  };
  return { light: of(FabChains.light.asym), heavy: of(FabChains.heavy.asym) };
}

/** The FASTA file: a record a chain, its name and the entity's description, the letters on one line. */
export function fabFasta(chains: ReturnType<typeof fabChains>): string {
  return (['light', 'heavy'] as const)
    .map((c) => `>${FabChains[c].name} ${chains[c].description}\n${chains[c].letters}\n`)
    .join('');
}

/** A numbered segment as SAbDab serves it: its letters and a `[[number, insertion], residue]` an IMGT position. */
export interface Segment {
  readonly sequence: string;
  readonly numbering: readonly (readonly [readonly [number, string], string])[];
}

/**
 * The segment's IMGT labels on its chain: the domain found in the chain's letters, a label its residue's (number and
 * insertion code), `null` outside the domain; a residue that is not the chain's there is refused.
 */
export function segmentLabels(segment: Segment, letters: string): (string | null)[] {
  const start = letters.indexOf(segment.sequence);
  if (start < 0) throw new Error('the segment is not in the chain');
  const out: (string | null)[] = Array.from(letters, () => null);
  let i = start;
  for (const [[n, ins], aa] of segment.numbering) {
    if (aa === '-') continue;
    if (letters[i] !== aa)
      throw new Error(`the segment's ${n}${ins.trim()} is ${aa}, the chain's ${i + 1} ${letters[i]}`);
    out[i++] = `${n}${ins.trim()}`;
  }
  return out;
}

/**
 * A numbering's `seq` file: a system a chain, a range a position with its `label`, missing where the position has none
 * (outside the variable domain: a source's numbering is a row a position).
 */
export function numberingSeq(chains: readonly { readonly name: string; readonly labels: Labels }[]): string {
  const parts: SeqData[] = chains.map((c) =>
    inlineSeq(
      { name: c.name, length: c.labels.length, alphabet: 'protein' },
      {
        ranges: {
          form: 'list',
          start: Int32Array.from(c.labels, (_, i) => i),
          end: Int32Array.from(c.labels, (_, i) => i + 1),
          fields: [{ name: 'label', column: categoryColumn(c.labels) }],
        },
      },
    ),
  );
  const data: SeqData = {
    coverage: {
      systems: parts.flatMap((p) => p.coverage.systems),
      shapes: parts.flatMap((p) => p.coverage.shapes),
      levels: [1],
    },
    pairs: parts.flatMap((p) => p.pairs),
  };
  return `${encodeSeq(data)}\n`;
}
