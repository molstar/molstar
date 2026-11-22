/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { columnValues } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { DefaultOpenLimits } from '../../data/protocol.js';
import { checkCoverage } from '../../data/seq-format.js';
import { GenbankRuntime, type GenbankSource } from '../genbank.js';
import { memoryAssets } from './bed-fixtures.js';

// The real INSDC record L09137.2 (pUC19c) as NCBI GenBank serves it (examples' `fixtures/plasmid/L09137.2.gb`, recorded
// as fetched): one cyclic record of 2,686 letters, the same letters as the recorded FASTA's, its one `source` feature,
// found by its LOCUS name and its accession; nothing refused but the header fields not read.

const fixture = (name: string) =>
  readFileSync(fileURLToPath(new URL(`../../../../examples/fixtures/plasmid/${name}`, import.meta.url)), 'utf8');
const signal = new AbortController().signal;

describe('the real L09137.2 GenBank record', () => {
  it('one cyclic record SYNPUC19CV (alias L09137.2) of 2,686 letters, one source item, found by name', async () => {
    const s = (await GenbankRuntime.open(
      { address: '/data/gb', format: 'genbank', uri: 'f' },
      { assets: memoryAssets({ f: fixture('L09137.2.gb') }), signal, limits: DefaultOpenLimits },
    )) as GenbankSource;
    expect(checkCoverage(s.coverage)).toBeUndefined();
    expect(s.coverage.systems).toEqual([
      { name: 'SYNPUC19CV', length: 2686, alphabet: 'dna', cyclic: true, aliases: ['L09137.2'] },
    ]);
    expect(s.records[0]).toMatchObject({
      definition: 'Cloning vector pUC19c, complete sequence.',
      organism: 'Cloning vector pUC19c',
      features: 1,
    });
    expect((s.issues ?? []).map((i) => i.code)).toEqual(['seq.source.genbank.header']);
    const a = await s.read(
      { system: 'SYNPUC19CV', start: 0, end: 2686, level: 1, shapes: ['letters', 'ranges'] },
      signal,
    );
    if (!a.ok || a.tables.ranges?.form !== 'list') throw new Error('not read');
    // the letters as read (lowercase), the recorded FASTA's (converted to the `seq` file) whatever their case
    const seq = JSON.parse(fixture('L09137.2.seq')).params.pairs[0].answer.tables.letters.letters as string;
    expect((a.tables.letters?.letters as string).toUpperCase()).toBe(seq);
    expect(columnValues(a.tables.ranges.id!)).toEqual(['source@1-2686']);
    const hits = async (q: string) => (await s.findNames!(q, { max: 5 }, signal)).map((h) => [h.field, h.system]);
    expect(await hits('synpuc19cv')).toEqual([['locus', 'SYNPUC19CV']]);
    expect(await hits('L09137.2')).toEqual([['accession', 'SYNPUC19CV']]);
  });
});
