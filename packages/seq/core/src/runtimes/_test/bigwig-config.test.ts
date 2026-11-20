/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { DefaultBigWigLimits } from '../bigwig.js';

// The bigWig runtime's Configuration table (_spec/README.md "Data", bigWig) against `DefaultBigWigLimits`: a row for
// every bound, and no row for one it lacks.

describe('the bigWig runtime’s Configuration table', () => {
  it('names every item of DefaultBigWigLimits, and none it lacks', () => {
    const md = readFileSync(fileURLToPath(new URL('../../../_spec/README.md', import.meta.url)), 'utf8');
    const section = md.slice(md.indexOf('\n### bigWig\n'), md.indexOf('\n### Archives\n'));
    const rows = [...section.matchAll(/^ *\| `(\w+)` +\|/gm)].map((m) => m[1]);
    expect(rows.sort()).toEqual(Object.keys(DefaultBigWigLimits).sort());
  });
});
