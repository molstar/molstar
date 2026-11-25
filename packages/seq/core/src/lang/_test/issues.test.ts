/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { BoundIssues, LangIssues } from '../issues.js';

describe('seq.lang: issue codes', () => {
  it('every code a source uses is cataloged, and every cataloged code is used', () => {
    const known = new Set([...Object.keys(LangIssues.codes), ...Object.keys(BoundIssues.codes)]);
    for (const info of Object.values(LangIssues.codes)) expect(info.severity).toMatch(/^(error|warning|info)$/);
    const dir = fileURLToPath(new URL('../', import.meta.url));
    const used = new Set<string>();
    for (const name of readdirSync(dir)) {
      if (!name.endsWith('.ts') || name === 'issues.ts') continue;
      for (const m of readFileSync(join(dir, name), 'utf8').matchAll(/'(seq\.lang\.[a-z_.]*[a-z_])'/g)) used.add(m[1]);
    }
    expect([...used].filter((c) => !known.has(c))).toEqual([]);
    expect([...known].filter((c) => !used.has(c))).toEqual([]);
  });
});
