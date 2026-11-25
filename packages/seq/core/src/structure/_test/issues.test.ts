/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { StructureSyncIssues } from '../issues.js';

describe('seq.structure: issue codes', () => {
  it('every code a source uses is cataloged, every cataloged code is used, and the spec lists each', () => {
    const known = new Set(Object.keys(StructureSyncIssues.codes));
    const dir = fileURLToPath(new URL('../', import.meta.url));
    const used = new Set<string>();
    for (const name of readdirSync(dir)) {
      if (!name.endsWith('.ts') || name === 'issues.ts') continue;
      for (const m of readFileSync(join(dir, name), 'utf8').matchAll(/'(seq\.structure\.[a-z_.]*[a-z_])'/g))
        used.add(m[1]);
    }
    expect([...used].filter((c) => !known.has(c))).toEqual([]);
    expect([...known].filter((c) => !used.has(c))).toEqual([]);
    const spec = readFileSync(fileURLToPath(new URL('../../../_spec/README.md', import.meta.url)), 'utf8');
    const listed = /`seq\.structure\.\*` \(`StructureSyncIssues`\): ([^]*?)\n\n/.exec(spec)![1];
    for (const c of known) expect(listed).toContain(`\`${c.slice('seq.structure.'.length)}\``);
  });
});
