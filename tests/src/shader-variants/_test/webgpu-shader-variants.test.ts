/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import {
  buildSources,
  compileSource,
  MVSDocuments,
  prepareAllKernels,
  recordCompiles,
  SourceNames,
  walkVariants,
} from '../matrix.js';
import { createVariantPlugin, listMVSDocuments, readRepositoryFile, WebGPU } from './node.js';

/** Programs Node has no means for: SMAA reads its lookup tables as `HTMLImageElement`s (checked in the browser) */
const NodeUnsupported = new Set(['smaa-edges', 'smaa-weights', 'smaa-blend']);

const takeErrors = (graphics: unknown) => (graphics as { takeErrors?: () => string[] }).takeErrors?.() ?? [];

// Minutes per source: on demand only (`SHADER_VARIANTS=1 pnpm native:run -- pnpm test:unit tests/src/shader-variants`)
const OnDemand = process.env.SHADER_VARIANTS === '1';

describe.skipIf(!OnDemand)('shader variants', () => {
  it('read every MolViewSpec document of the repository', () => {
    const documents = [
      ...listMVSDocuments('data/examples/mvs/'),
      ...listMVSDocuments('packages/mvs/runtime/src/_test/fixtures/'),
    ];
    expect([...MVSDocuments].sort()).toEqual(documents.sort());
  });
});

// The variants scenes reach, compiled on WebGPU: each source alone (`-t '<source>'` runs one) without an error or a
// warning from Dawn (Tint), the kernels, and all of them for the coverage of every program and kernel
describe.skipIf(!OnDemand || !WebGPU)('shader variants on WebGPU', () => {
  it.each(SourceNames)(
    '%s compile cleanly',
    async (name) => {
      const { plugin, dispose } = await createVariantPlugin({ webgpu: WebGPU! }, 'webgpu');
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
    },
    600_000,
  );

  it('kernels compile cleanly', async () => {
    const { plugin, dispose } = await createVariantPlugin({ webgpu: WebGPU! }, 'webgpu');
    try {
      const { graphics } = plugin.canvas3d!;
      const record = recordCompiles(graphics);
      await prepareAllKernels(graphics);
      expect(await record.problems()).toEqual([]);
      expect(takeErrors(graphics)).toEqual([]);
      record.stop();
    } finally {
      dispose();
    }
  }, 120_000);

  it('cover every program and kernel', async () => {
    const { plugin, dispose } = await createVariantPlugin({ webgpu: WebGPU! }, 'webgpu');
    try {
      const { graphics } = plugin.canvas3d!;
      const record = recordCompiles(graphics);
      const walk = await walkVariants(plugin, readRepositoryFile, record.count);
      await prepareAllKernels(graphics);

      // programs first compiled by preparing frames (passes precompile does not plan), and by drawing them (passes
      // preparing does not plan)
      const precompiled = record.compiledIn(0, walk.precompiled);
      const prepared = record.compiledIn(0, walk.prepared);
      const byPrepare = [...record.compiledIn(walk.precompiled, walk.prepared)].filter((n) => !precompiled.has(n));
      const byFrames = [...record.compiledIn(walk.prepared, walk.drawn)].filter((n) => !prepared.has(n));
      console.log(
        `shader variants on WebGPU\n  ${walk.summary.join('\n  ')}\n  first by preparing: ${byPrepare.join(', ')}; ` +
          `first by drawing: ${byFrames.join(', ')}`,
      );
      expect(record.uncovered().filter((name) => !NodeUnsupported.has(name))).toEqual([]);
      record.stop();
    } finally {
      dispose();
    }
  }, 600_000);
});
