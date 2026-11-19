/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { execFileSync } from 'node:child_process';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';

const script = join(dirname(fileURLToPath(import.meta.url)), '../sizes.mjs');

describe('the size script', () => {
  // the smallest entry: a bundle of the sources, offline
  it('--files lists the largest input files of each entry, as its usage line says', () => {
    const out = execFileSync('node', [script, '--only', 'spec-read.ts', '--files', '--top', '3'], {
      encoding: 'utf8',
      timeout: 60_000,
    });
    expect(out).toMatch(/reading a spec[\s\S]*\n {2}files:\n( {4}\d+\.\d kB\s+\S+\n){3}/);
  }, 60_000);
});
