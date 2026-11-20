/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { existsSync, readFileSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { syntheticSignal, writeBigWig } from './bigwig-writer.js';

// p53's conservation signal over the gene's region, a bigWig the examples commit (`fixtures/p53/`) and read by region:
// generated here, where the test-only writer is, from a seed (never recorded from a service). The examples' own test
// reads it through the runtime against the region's system and length. `SEQ_RECORDS=update` writes it again.

/** The examples' gene region (`p53/region.ts`: `RegionSystem`, `RegionLength`). */
const System = 'chr17-tp53';
const Length = 6000;
const File = fileURLToPath(new URL('../../../../examples/fixtures/p53/gene-region-conservation.bw', import.meta.url));

/** Conservation a base in [0, 1), in runs; the region's last 200 bases without data (past the alignment's end). */
function p53Conservation() {
  return writeBigWig([{ name: System, values: syntheticSignal(53, Length, [[Length - 200, Length]]) }], {
    zooms: [16, 256],
    itemsPerSection: 256,
    blockSize: 16,
    // deflate's bytes differ by zlib build (Node versions): stored as they are, the file is the same everywhere
    compress: false,
  });
}

describe('p53’s conservation bigWig', () => {
  it('the committed file is what its generator writes', async () => {
    const bytes = await p53Conservation();
    if (process.env.SEQ_RECORDS === 'update') writeFileSync(File, bytes);
    expect(existsSync(File), `no ${File}: run with SEQ_RECORDS=update`).toBe(true);
    expect(new Uint8Array(readFileSync(File))).toEqual(bytes);
  });
});
