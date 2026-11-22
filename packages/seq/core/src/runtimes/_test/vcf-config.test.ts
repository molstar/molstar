/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { DefaultVcfOptions } from '../vcf.js';

// The VCF runtime's Configuration table (_spec/README.md "Data", VCF) against `DefaultVcfOptions`: a row for every
// bound and for `tabix` (io's other bounds, optional), and no row for one it lacks.

describe('the VCF runtime’s Configuration table', () => {
  it('names every item of DefaultVcfOptions and `tabix`, and none it lacks', () => {
    const md = readFileSync(fileURLToPath(new URL('../../../_spec/README.md', import.meta.url)), 'utf8');
    const start = md.indexOf('\n### VCF\n');
    const section = md.slice(start, md.indexOf('\n### ', start + 1));
    const rows = [...section.matchAll(/^ *\| `(\w+)` +\|/gm)].map((m) => m[1]);
    expect(rows.sort()).toEqual([...Object.keys(DefaultVcfOptions), 'tabix'].sort());
  });
});
