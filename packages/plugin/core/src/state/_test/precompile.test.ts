/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { createVolumeRepresentationParams } from '@molstar/plugin/state/helpers/volume-representation-params';
import { getTreePrecompileList } from '@molstar/plugin/state/precompile';
import { VolumeRepresentation3D } from '@molstar/plugin/state/transforms/volume/representation';

describe('precompile lists from state trees', () => {
  it('take an undefined param as its default (MVS writes `clip: undefined` without clip nodes)', async () => {
    const plugin = new PluginContext(DefaultPluginSpec());
    await plugin.init();
    try {
      const params = createVolumeRepresentationParams(plugin, undefined, { type: 'isosurface' });
      const derive = (typeParams: object) => {
        const builder = plugin.build();
        builder.toRoot().apply(VolumeRepresentation3D, { ...params, type: { name: 'isosurface', params: typeParams } });
        return getTreePrecompileList(plugin, builder.getTree()).objects;
      };
      const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
      const withDefaults = derive({ ...params.type.params });
      const withUndefined = derive({ ...params.type.params, clip: undefined });
      expect(warn).not.toHaveBeenCalled();
      expect(withUndefined.length).toBeGreaterThan(0);
      expect(withUndefined).toEqual(withDefaults);
    } finally {
      plugin.dispose();
    }
  });
});
