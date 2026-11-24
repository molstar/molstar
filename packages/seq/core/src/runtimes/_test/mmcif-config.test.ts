/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import type { MmcifOptions } from '../mmcif.js';

// The mmCIF runtime's Configuration table (_spec/README.md "Data", mmCIF) against `MmcifOptions`: a row for every
// option, and no row for one it lacks.

const Options: Record<keyof MmcifOptions, true> = { assembly: true, water: true, name: true, limits: true };

describe('the mmCIF runtime’s Configuration table', () => {
  it('names every option of MmcifOptions, and none it lacks', () => {
    const md = readFileSync(fileURLToPath(new URL('../../../_spec/README.md', import.meta.url)), 'utf8');
    const start = md.indexOf('\n### mmCIF\n');
    const section = md.slice(start, md.indexOf('\n### ', start + 1));
    const rows = [...section.matchAll(/^ *\| `(\w+)` +\|/gm)].map((m) => m[1]);
    expect(rows.sort()).toEqual(Object.keys(Options).sort());
  });
});
