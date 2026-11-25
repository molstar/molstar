/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { expectOk } from '../../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../../default-spec.js';
import { createPrismViewerHeadless } from '../../../headless.js';
import type { PrismViewer } from '../../../viewer.js';

// `canvas.backend` (the plugin's `setGraphicsBackends`) on a stand-in canvas: there is no WebGL or WebGPU here.

let open: PrismViewer | undefined;
afterEach(async () => {
  await open?.dispose();
  open = undefined;
});

/** A viewer whose plugin was made from `ids` (a WebGPU set lists WebGL's too), drawn with the first; `creatable` can be made. */
async function withBackends(ids: readonly string[], creatable: readonly string[] = ids) {
  const v = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec(), modules: ['canvas'] });
  open = v;
  const capabilities = { backend: ids[0] };
  const switched: string[] = [];
  Object.defineProperty(v.plugin, 'canvas3d', {
    value: { graphics: { capabilities }, dispose() {} },
    configurable: true,
  });
  Object.defineProperty(v.plugin, 'graphicsBackends', {
    value: [...ids, ...ids].map((id) => ({ id })),
    configurable: true,
  });
  Object.defineProperty(v.plugin, 'setGraphicsBackends', {
    value: async (backends: readonly { id: string }[]) => {
      const next = backends.find((b) => creatable.includes(b.id));
      if (!next) return false;
      switched.push(next.id);
      capabilities.backend = next.id;
      return true;
    },
    configurable: true,
  });
  return { v, switched };
}

describe('canvas.backend', () => {
  it('switches to a backend the viewer has, no history step; already-so', async () => {
    const { v, switched } = await withBackends(['webgpu', 'webgl2', 'webgl']);
    expect(v.api.describe('canvas.backend')?.available).toBe(true);
    const r = await v.api.canvas.backend!({ backend: 'webgl2' });
    expectOk(r);
    expect(r.ok && r.did).toBe('Switched to WebGL 2.');
    expect(switched).toEqual(['webgl2']);
    expect(await v.api.canvas.backend!({ backend: 'webgl2' })).toMatchObject({ ok: true, unchanged: true });
    expect(switched).toEqual(['webgl2']);
    expect(v.managers.history.steps.value).toHaveLength(0);
  });

  it('says which backends there are, and when one cannot be created', async () => {
    const { v } = await withBackends(['webgl2', 'webgl'], ['webgl2']);
    expect(await v.api.canvas.backend!({ backend: 'webgpu' })).toMatchObject({ ok: false, code: 'invalid' });
    expect(await v.api.canvas.backend!({ backend: 'webgl' })).toMatchObject({ ok: false, code: 'failed' });
  });

  it('is not offered with one backend', async () => {
    const { v } = await withBackends(['webgl2']);
    expect(v.api.describe('canvas.backend')?.available).toMatch(/one graphics backend/);
  });
});
