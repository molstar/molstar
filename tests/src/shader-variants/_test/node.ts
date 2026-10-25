/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as fs from 'node:fs';
import { readFile } from 'node:fs/promises';
import { createRequire } from 'node:module';
import { pathToFileURL } from 'node:url';
import { setFSModule } from '@molstar/core/util/data-source';
import { clearFontAtlasCache, setCanvasModule } from '@molstar/graphics/geo/geometry/text/font-atlas';
import { MolViewSpec } from '@molstar/mvs/behavior';
import { MVSRuntimeRegistry } from '@molstar/mvs/registry';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import type { ExternalModules } from '@molstar/plugin-headless/screenshot';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { PluginSpec } from '@molstar/plugin/spec';
import type { ReadFile } from '../matrix.js';

// The Node side of the shader variants (`.v6/plans/gl/shader-analysis.md`): the repository's files, a headless plugin,
// the optional native modules. Test helpers, not a test file itself.

setFSModule(fs);

const root = new URL('../../../../', import.meta.url);

/** A file of the repository by its path from the root */
export const readRepositoryFile: ReadFile = (path) => readFile(new URL(path, root), 'utf8');

/** The `.mvsj` files of a directory of the repository, by path from the root */
export function listMVSDocuments(dir: string) {
  return fs
    .readdirSync(new URL(dir, root))
    .filter((f) => f.endsWith('.mvsj'))
    .map((f) => `${dir}${f}`);
}

const require = createRequire(import.meta.url);

/** The optional `webgpu` package (Dawn), when it loads and has an adapter */
export const WebGPU = await (async () => {
  try {
    let module: NonNullable<ExternalModules['webgpu']>;
    try {
      module = (await import('webgpu' as string)) as NonNullable<ExternalModules['webgpu']>;
    } catch {
      module = (await import(pathToFileURL(require.resolve('webgpu')).href)) as NonNullable<ExternalModules['webgpu']>;
    }
    return (await module.create([]).requestAdapter()) ? module : undefined;
  } catch {
    return undefined;
  }
})();

/** The optional native `gl` package (WebGL1), when it loads and makes a context */
export const GL = (() => {
  try {
    const module = require('gl') as typeof import('gl');
    const context = module(1, 1);
    if (!context) return undefined;
    context.getExtension('STACKGL_destroy_context')?.destroy();
    return module;
  } catch {
    return undefined;
  }
})();

/**
 * A headless plugin with the MolViewSpec runtime (for its documents) on `backend`, its canvas large enough for every
 * level of the blurs and pyramids; with a font that draws nothing (labels draw their glyphs with a 2D canvas, which
 * Node has not, and only the programs matter here). `dispose` releases both
 */
export async function createVariantPlugin(modules: ExternalModules, backend: 'webgl' | 'webgpu') {
  const context = {
    measureText: (text: string) => ({ width: 10 * text.length, actualBoundingBoxLeft: 0, actualBoundingBoxRight: 10 }),
    clearRect() {},
    fillText() {},
    getImageData: (_x: number, _y: number, w: number, h: number) => ({ data: new Uint8ClampedArray(w * h * 4) }),
  };
  setCanvasModule({ createCanvas: () => ({ getContext: () => context }) });
  const spec = DefaultPluginSpec();
  const plugin = await HeadlessPluginContext.create(
    modules,
    {
      ...spec,
      registry: [...(spec.registry ?? []), MVSRuntimeRegistry],
      behaviors: [...spec.behaviors, PluginSpec.Behavior(MolViewSpec)],
    },
    { width: 256, height: 256 },
    { backend },
  );
  await plugin.init();
  return {
    plugin,
    dispose: () => {
      plugin.dispose();
      setCanvasModule(undefined);
      clearFontAtlasCache();
    },
  };
}
