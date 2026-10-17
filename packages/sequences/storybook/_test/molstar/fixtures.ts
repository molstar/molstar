/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Mol* demo's fixtures for the tests of the demo's support code (src/_stories/support): the structures are the
// Mol* adapter's recorded files (Prism's, packages/prism/viewer/src/sequence/molstar/_test/fixtures, read as data), the
// p53 data the demo's (../fixtures/p53), P04637's sequence Seq*'s shared fixture (core's p53.json). No network.

import { readFileSync } from 'node:fs';
import { createRegistries, KernelRegistry, type Loaders, MvsRegistry, register } from '@molstar/seq-core/format';
import { MolstarSvsEntry, structureUrl } from '@molstar/prism-viewer/sequence/molstar';
import { P53_STRUCTURES, p53Loaders } from '../../src/_stories/support/seqstar-p53-demo.js';

const STRUCTURES = '../../../../prism/viewer/src/sequence/molstar/_test/fixtures/structures/';
export const bytesOf = (file: string) => new Uint8Array(readFileSync(new URL(`${STRUCTURES}${file}`, import.meta.url)));
export const P53: { sequence: string } = JSON.parse(
  readFileSync(new URL('../../../core/_test/fixtures/p53.json', import.meta.url), 'utf8'),
);

const fasta = (header: string, seq: string) => `>${header}\n${seq.replace(/(.{60})/g, '$1\n')}\n`;

/** The demo's loaders (seqstar-p53-demo.ts) over the committed fixtures. */
export function demoLoaders(): Loaders {
  const structures = new Map<string, Uint8Array>();
  for (const s of P53_STRUCTURES) {
    const at = structureUrl({ source: s.source });
    if (at)
      structures.set(
        at.url,
        bytesOf(`${s.source.id.toLowerCase() === 'af-p04637-f1' ? 'AF-P04637-F1-model_v6' : s.source.id}.bcif`),
      );
  }
  return p53Loaders({
    structures,
    fasta: fasta('sp|P04637|P53_HUMAN Cellular tumor antigen p53', P53.sequence),
    features: readFileSync(new URL('../fixtures/p53/P04637.json', import.meta.url), 'utf8'),
    variants: readFileSync(new URL('../fixtures/p53/variation-P04637.json', import.meta.url), 'utf8'),
    ptms: readFileSync(new URL('../fixtures/p53/P04637-ptms.tsv', import.meta.url), 'utf8'),
  });
}

export function registries() {
  const reg = createRegistries();
  register(reg, [KernelRegistry, MvsRegistry, MolstarSvsEntry]);
  return reg;
}
