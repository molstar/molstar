/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { Asset, DefaultAssetManager } from '@molstar/core/util/assets';
import { PluginContext } from '@molstar/plugin/context';

function fileAsset() {
  return Asset.File(new File(['abc'], 'a.txt'));
}

describe('PluginSpec.assetManager', () => {
  it('uses a new DefaultAssetManager when undefined, and disposes it', () => {
    const plugin = new PluginContext({ behaviors: [] });
    const manager = plugin.managers.asset;
    expect(manager).toBeInstanceOf(DefaultAssetManager);
    expect(new PluginContext({ behaviors: [] }).managers.asset).not.toBe(manager);

    const asset = fileAsset();
    manager.set(asset, asset.file!);
    plugin.dispose();
    expect(manager.has(asset)).toBe(false);
  });

  it('calls the factory once per plugin and uses what it returns', async () => {
    const manager = new DefaultAssetManager();
    const factory = vi.fn(() => ({ manager }));
    const plugin = new PluginContext({ behaviors: [], assetManager: factory });
    expect(factory).toHaveBeenCalledTimes(1);
    expect(plugin.managers.asset).toBe(manager);

    const asset = fileAsset();
    const wrapper = await plugin.runTask(plugin.managers.asset.resolve(asset, 'string'));
    expect(wrapper.data).toBe('abc');
    expect(manager.has(asset)).toBe(true);
    wrapper.dispose();
    expect(manager.has(asset)).toBe(false);
  });

  it('disposes the manager with the plugin unless disposeWithPlugin is false', () => {
    const owned = new DefaultAssetManager();
    const shared = new DefaultAssetManager();
    const ownedDispose = vi.spyOn(owned, 'dispose');
    const sharedDispose = vi.spyOn(shared, 'dispose');

    const a = new PluginContext({ behaviors: [], assetManager: () => ({ manager: owned }) });
    const b = new PluginContext({ behaviors: [], assetManager: () => ({ manager: shared, disposeWithPlugin: false }) });
    const c = new PluginContext({ behaviors: [], assetManager: () => ({ manager: shared, disposeWithPlugin: false }) });
    expect(b.managers.asset).toBe(c.managers.asset);

    a.dispose();
    b.dispose();
    c.dispose();
    expect(ownedDispose).toHaveBeenCalledTimes(1);
    expect(sharedDispose).not.toHaveBeenCalled();
  });
});
