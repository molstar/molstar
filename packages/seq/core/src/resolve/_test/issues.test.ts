/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { GenomeIssues } from '../genome-issues.js';
import { ResolveIssues } from '../issues.js';

describe('seq.resolve: issue codes', () => {
  it('every code a source uses is cataloged, and every cataloged code is used', () => {
    const known = new Set([...Object.keys(ResolveIssues.codes), ...Object.keys(GenomeIssues.codes)]);
    const dir = fileURLToPath(new URL('../', import.meta.url));
    const used = new Set<string>();
    for (const name of readdirSync(dir)) {
      if (!name.endsWith('.ts') || name === 'issues.ts' || name === 'genome-issues.ts') continue;
      for (const m of readFileSync(join(dir, name), 'utf8').matchAll(/'(seq\.resolve\.[a-z_.]*[a-z_])'/g))
        used.add(m[1]);
    }
    expect([...used].filter((c) => !known.has(c))).toEqual([]);
    expect([...known].filter((c) => !used.has(c))).toEqual([]);
  });
});
