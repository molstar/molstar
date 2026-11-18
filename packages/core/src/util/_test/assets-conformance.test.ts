/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { Task } from '@molstar/core/task';
import { type Asset, type AssetManager, DefaultAssetManager } from '../assets.js';
import type { DataType } from '../data-source.js';
import { AssetManagerCases, checkAssetManager } from '../assets-conformance.js';

/** What an app writes: the default manager with credentials of its own when a read names none. */
class CredentialedManager extends DefaultAssetManager {
  override resolve<T extends DataType>(
    asset: Asset,
    type: T,
    options: boolean | AssetManager.ResolveOptions = true,
  ): Task<Asset.Wrapper<T>> {
    const given = typeof options === 'boolean' ? { store: options } : options;
    return super.resolve(asset, type, { credentials: 'include', ...given });
  }
}

/** A manager that drops every network option: the suite must tell. */
class IgnoringManager extends DefaultAssetManager {
  override resolve<T extends DataType>(
    asset: Asset,
    type: T,
    options: boolean | AssetManager.ResolveOptions = true,
  ): Task<Asset.Wrapper<T>> {
    const given = typeof options === 'boolean' ? { store: options } : options;
    const { credentials, redirect, signal, maxBytes, ...rest } = given;
    return super.resolve(asset, type, rest);
  }
}

describe('checkAssetManager', () => {
  it('has cases', () => {
    expect(AssetManagerCases.length).toBeGreaterThan(5);
  });

  it('DefaultAssetManager conforms', async () => {
    expect(await checkAssetManager(() => new DefaultAssetManager())).toEqual([]);
  });

  it('a derived manager that adds credentials conforms', async () => {
    expect(await checkAssetManager(() => new CredentialedManager())).toEqual([]);
  });

  it('a manager that ignores the four options is caught on each', async () => {
    const failures = await checkAssetManager(() => new IgnoringManager());
    for (const name of ['credentials are honored', 'redirects are honored', 'a signal cancels', 'maxBytes bounds']) {
      expect(
        failures.some((f) => f.startsWith(name)),
        name,
      ).toBe(true);
    }
    expect(failures.some((f) => f.startsWith('references are counted'))).toBe(false);
  });
});
