/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { buildSources, compileSource, recordCompiles, SourceNames, walkVariants } from '../matrix.js';
import { createVariantPlugin, GL, readRepositoryFile } from './node.js';

/** Programs the native `gl` package can't run: SMAA reads its lookup tables as `HTMLImageElement`s */
const NodeUnsupported = new Set(['smaa-edges', 'smaa-weights', 'smaa-blend']);

const takeErrors = (graphics: unknown) => (graphics as { takeErrors?: () => string[] }).takeErrors?.() ?? [];

// The variants scenes reach, compiled and linked as WebGL1 GLSL by the native `gl` package: each source alone
// (`-t '<source>'` runs one) without a compile or link failure, and all of them for the coverage of every program.
// Minutes per source: on demand only (`SHADER_VARIANTS=1 pnpm native:run -- pnpm test:unit tests/src/shader-variants`)
describe.skipIf(process.env.SHADER_VARIANTS !== '1' || !GL)('shader variants on WebGL', () => {
  for (const name of SourceNames) {
    it(`${name} compile and link`, async () => {
      const { plugin, dispose } = await createVariantPlugin({ gl: GL! }, 'webgl');
      try {
        const { graphics } = plugin.canvas3d!;
        const record = recordCompiles(graphics);
        const sources = await buildSources(plugin, readRepositoryFile);
        console.log(await compileSource(plugin, sources.find((s) => s.name === name)!, record.count));
        expect(await record.problems()).toEqual([]);
        expect(takeErrors(graphics)).toEqual([]);
        record.stop();
      } finally {
        dispose();
      }
    }, 900_000);
  }

  it('cover every program', async () => {
    const { plugin, dispose } = await createVariantPlugin({ gl: GL! }, 'webgl');
    try {
      const record = recordCompiles(plugin.canvas3d!.graphics);
      const walk = await walkVariants(plugin, readRepositoryFile, record.count);
      console.log(`shader variants on WebGL\n  ${walk.summary.join('\n  ')}`);
      expect(record.uncovered().filter((name) => !NodeUnsupported.has(name))).toEqual([]);
      record.stop();
    } finally {
      dispose();
    }
  }, 900_000);
});
