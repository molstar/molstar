/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { StateTransformer } from '@molstar/core/state';
import { CreateVolumeStreamingBehavior } from '../transformers.js';

const params = (name: string, channels: object) => ({ entry: { name, params: { channels, view: { name: 'auto' } } } });

const behavior = (info: object, entries: string[]) => {
  const update = vi.fn(async () => true);
  const b = {
    data: {
      data: info,
      infoMap: new Map(entries.map((e) => [e, { dataId: e }])),
      update,
      getDescription: () => '',
    },
  };
  return { b, update };
};

const runUpdate = (a: object, b: object, oldParams: object, newParams: object) =>
  (CreateVolumeStreamingBehavior.definition.update as any)({ a, b, oldParams, newParams }).run();

describe('the volume streaming behavior transform', () => {
  it('is recreated when its info is new (another state streams other maps)', async () => {
    // An x-ray entry's behavior, then a state that streams an EM map: the old behavior has no entry of that name.
    const { b, update } = behavior({}, ['1iep']);
    const newInfo = { data: {} };
    const result = await runUpdate(
      newInfo,
      b,
      params('1iep', { '2fo-fc': {} }),
      params('EMD-21375', { em: { isoValue: {} } }),
    );
    expect(result).toBe(StateTransformer.UpdateResult.Recreate);
    expect(update).not.toHaveBeenCalled();
  });

  it('is updated in place for the same info', async () => {
    const info = {};
    const { b, update } = behavior(info, ['EMD-1', 'EMD-2']);
    const result = await runUpdate({ data: info }, b, params('EMD-1', { em: {} }), params('EMD-2', { em: {} }));
    expect(result).toBe(StateTransformer.UpdateResult.Updated);
    expect(update).toHaveBeenCalledOnce();
  });
});
