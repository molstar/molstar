/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { PluginConfig } from '@molstar/plugin/config';
import { resolveOptions } from '../../options.js';
import type { PartialOptions } from '../../spec.js';
import { optionsConfig } from '../compose-plugin.js';

function backends(options: PartialOptions) {
  const item = optionsConfig(resolveOptions(options)).find(([i]) => i === PluginConfig.General.GraphicsBackends);
  return (item?.[1] as { id: string }[] | undefined)?.map((b) => b.id);
}

describe('optionsConfig: graphics backends', () => {
  it('lists WebGPU first by default (auto), then WebGL', () => {
    expect(backends({})).toEqual(['webgpu', 'webgl2', 'webgl']);
  });

  it('lists WebGL2, then WebGL1 for webgl, one version for webgl2 and webgl1', () => {
    expect(backends({ renderer: { graphicsBackend: 'webgl' } })).toEqual(['webgl2', 'webgl']);
    expect(backends({ renderer: { graphicsBackend: 'webgl2' } })).toEqual(['webgl2']);
    expect(backends({ renderer: { graphicsBackend: 'webgl1' } })).toEqual(['webgl']);
  });
});
