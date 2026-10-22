/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import path from 'node:path';
import { fileURLToPath } from 'node:url';
import * as esbuild from 'esbuild';
import { describe, expect, it } from 'vitest';
import { startBytes } from '../start-closure.mjs';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../../..');

describe('ESM start closure (perf audit 2, 12/B2)', () => {
  const out = (bytes, ...imports) => ({ bytes, imports });
  const meta = {
    outputs: {
      'build/index.js': out(10, { path: 'build/chunks/a.js', kind: 'import-statement' }),
      'build/chunks/a.js': out(
        100,
        { path: 'build/chunks/b.js', kind: 'import-statement' },
        { path: 'build/chunks/lazy.js', kind: 'dynamic-import' },
      ),
      'build/chunks/b.js': out(
        1000,
        { path: 'build/chunks/a.js', kind: 'import-statement' },
        { path: 'crypto', kind: 'import-statement', external: true },
      ),
      'build/chunks/lazy.js': out(10000, { path: 'build/chunks/b.js', kind: 'import-statement' }),
    },
  };

  it('counts the entry and its static chunks once, not lazy ones', () => {
    expect(startBytes(meta, 'build/index.js')).toBe(1110);
    expect(startBytes(meta, 'build/chunks/lazy.js')).toBe(11100);
  });

  it('throws on an entry the metafile does not have, instead of counting 0 bytes', () => {
    expect(() => startBytes(meta, 'build/other.js')).toThrow(/no output/);
  });

  it("keeps Seq* out of the AFDB PAE reader's imports (12/P9)", async () => {
    // the reader is in Prism's start code (model-archive's PAE operations); the Seq* barrel would add ~1.1 MB to it
    const r = await esbuild.build({
      absWorkingDir: root,
      entryPoints: ['packages/extensions/model-archive/src/pae/afdb.ts'],
      bundle: true,
      write: false,
      metafile: true,
      format: 'esm',
      conditions: ['molstar-src', 'import', 'default'],
      logLevel: 'silent',
    });
    const inputs = Object.keys(r.metafile.inputs);
    expect(inputs.filter((f) => f.startsWith('packages/sequences/') || f.includes('sequence/molstar/index'))).toEqual(
      [],
    );
    expect(inputs).toContain('packages/prism/viewer/src/sequence/molstar/pae-matrix.ts');
  });
});
