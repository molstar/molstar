/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { expect, it, vi } from 'vitest';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { resolveVideoParams, type VideoRendererBackend } from '../backend.js';

const config = { width: 64, height: 64, fps: 30 };
function fixture() {
  const params = {
    container: PD.Select('mp4', [['mp4', 'MP4']]),
    codec: PD.MappedStatic('avc', {
      avc: PD.Group({
        hardwareAcceleration: PD.Select('no-preference', [
          ['no-preference', 'Automatic'],
          ['prefer-hardware', 'Hardware'],
        ]),
        rate: PD.MappedStatic('high', {
          high: PD.EmptyGroup(),
          bitrate: PD.Group({
            bitrate: PD.Numeric(8),
            mode: PD.Select('variable', [
              ['variable', 'Variable'],
              ['constant', 'Constant'],
            ]),
          }),
        }),
      }),
    }),
  };
  const validate = vi.fn(async () => undefined);
  const backend: VideoRendererBackend = {
    name: 'Test',
    dimensionAlignment: 1,
    getParams: async () => params,
    validate,
    create: async () => {
      throw new Error('Not used');
    },
  };
  return { backend, validate };
}

it.each([
  [{ hardwareAcceleration: 'prefer-hadrware' }, 'hardwareAcceleration'],
  [{ rate: { name: 'typo', params: {} } }, 'rate'],
  [{ rate: { name: 'bitrate', params: { mode: 'typo' } } }, 'mode'],
  [{ rate: { name: 'bitrate', params: { unknown: 8 } } }, 'unknown'],
])('rejects explicit invalid nested choices instead of substituting defaults', async (params, path) => {
  const f = fixture();
  await expect(resolveVideoParams(f.backend, config, { codec: { name: 'avc', params } })).rejects.toThrow(path);
  expect(f.validate).not.toHaveBeenCalled();
});

it('fills omitted defaults and allows shared FPS alongside backend settings', async () => {
  const f = fixture();
  const values = await resolveVideoParams(f.backend, config, {
    fps: 30,
    codec: { name: 'avc', params: { rate: { name: 'bitrate', params: { bitrate: 8.3 } } } },
  });
  expect(values.codec.params.rate.params).toEqual({ bitrate: 8.3, mode: 'variable' });
  expect(values.codec.params.hardwareAcceleration).toBe('no-preference');
  expect(values.fps).toBeUndefined();
});
