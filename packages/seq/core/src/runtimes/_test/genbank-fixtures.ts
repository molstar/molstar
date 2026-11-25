/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { PCG } from '@molstar/core/data/util/hash-functions';

// GenBank written by a generator (the runtime's tests read each): records of letters with features of every location
// form, and the examples' generated plasmid record (`fixtures/plasmid/pUC19c-rotated.gb`) as text.

/** The examples' generated pUC19c record (the rotated variant: bla a join over the origin, lacZα a complement). */
export function plasmidGenbank(): string {
  const path = fileURLToPath(new URL('../../../../examples/fixtures/plasmid/pUC19c-rotated.gb', import.meta.url));
  return readFileSync(path, 'utf8');
}

/** ORIGIN's lines of `letters` (lowercase, 60 a line in blocks of 10, numbered). */
export function originLines(letters: string): string[] {
  const out = ['ORIGIN'];
  for (let at = 0; at < letters.length; at += 60) {
    const blocks: string[] = [];
    for (let k = at; k < Math.min(at + 60, letters.length); k += 10) blocks.push(letters.slice(k, k + 10));
    out.push(`${String(at + 1).padStart(9)} ${blocks.join(' ')}`);
  }
  return out;
}

export interface GbFeature {
  readonly key: string;
  readonly location: string;
  readonly qualifiers?: readonly (readonly [string, string?])[];
}

/** A record's text: LOCUS, DEFINITION, ACCESSION and VERSION, SOURCE, the features, ORIGIN, `//`. */
export function record(
  name: string,
  letters: string,
  features: readonly GbFeature[],
  options: { circular?: boolean; molecule?: string; accession?: string; length?: number; protein?: boolean } = {},
): string {
  const unit = options.protein ? 'aa' : 'bp';
  const molecule = options.protein ? '' : (options.molecule ?? 'DNA');
  const lines = [
    `LOCUS       ${name.padEnd(16)} ${String(options.length ?? letters.length).padStart(11)} ${unit}    ${molecule.padEnd(7)} ${options.circular ? 'circular' : 'linear  '} SYN 01-JAN-2026`,
    `DEFINITION  A generated record ${name}.`,
    `ACCESSION   ${options.accession ?? `X${name}`}`,
    `VERSION     ${options.accession ?? `X${name}`}.1`,
    'SOURCE      synthetic construct',
    '  ORGANISM  synthetic construct',
    '            other sequences; artificial sequences.',
    'FEATURES             Location/Qualifiers',
  ];
  for (const f of features) {
    lines.push(`     ${f.key.padEnd(16)}${f.location}`);
    for (const [k, v] of f.qualifiers ?? []) lines.push(`                     /${k}${v === undefined ? '' : `=${v}`}`);
  }
  lines.push(...originLines(letters), '//');
  return `${lines.join('\n')}\n`;
}

/** Letters of `n` nucleotides from a seed. */
export function letters(n: number, seed = 3): string {
  const rng = new PCG(seed);
  let s = '';
  for (let i = 0; i < n; i++) s += 'acgt'[Math.floor(rng.float() * 4)];
  return s;
}

/** A generated record of `features` genes on `n` bases: genes, joined CDSs on both strands, some partial. */
export function generated(features: number, n = features * 120 + 1000, seed = 7): string {
  const rng = new PCG(seed);
  const list: GbFeature[] = [];
  for (let i = 0; i < features; i++) {
    const a = 1 + Math.floor(rng.float() * (n - 400));
    const b = a + 100 + Math.floor(rng.float() * 250);
    const loc =
      i % 3 === 0
        ? `${a}..${b}`
        : i % 3 === 1
          ? `complement(join(${a}..${a + 20},${a + 40}..${b}))`
          : `join(<${a}..${a + 30},${a + 60}..>${b})`;
    list.push({
      key: i % 2 ? 'CDS' : 'gene',
      location: loc,
      qualifiers: [
        ['gene', `"g${i}"`],
        ['locus_tag', `"T_${i}"`],
        ['product', `"protein ${i}"`],
      ],
    });
  }
  return record('gen', letters(n, seed), list);
}
