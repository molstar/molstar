/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PresetProvider } from './preset-provider.js';

/**
 * Reference-counted registry of presets, keyed by `id` and (when present) `alias`.
 *
 * - Registering the same object again increments its count.
 * - Registering a different object under an existing key throws.
 * - Unregistering decrements the count and removes the preset at zero; an unknown preset is a no-op.
 */
export class PresetRegistry<P extends PresetProvider<any, any, any, any, any>> {
  private _providers: P[] = [];
  private keys = new Map<string, P>();
  private counts = new Map<P, number>();

  constructor(private readonly name: string) {}

  /** Providers in registration order. */
  get providers(): ReadonlyArray<P> {
    return this._providers;
  }

  /** The message `register` would throw for this provider, or `undefined` if it would succeed. */
  findConflict(provider: P): string | undefined {
    const conflict = (key: string, kind: 'id' | 'alias') => {
      const existing = this.keys.get(key);
      if (existing && existing !== provider) {
        return `${this.name}: ${kind} '${key}' of preset '${provider.id}' is already registered by a different preset ('${existing.id}')`;
      }
    };
    const byId = conflict(provider.id, 'id');
    if (byId) return byId;
    if (provider.alias !== undefined) return conflict(provider.alias, 'alias');
  }

  register(provider: P) {
    const conflict = this.findConflict(provider);
    if (conflict) throw new Error(conflict);

    const count = this.counts.get(provider);
    if (count !== undefined) {
      this.counts.set(provider, count + 1);
      return;
    }
    this.counts.set(provider, 1);
    this._providers.push(provider);
    this.keys.set(provider.id, provider);
    if (provider.alias !== undefined) this.keys.set(provider.alias, provider);
  }

  /** Accepts a provider or a provider id. */
  unregister(providerOrId: P | string) {
    const provider = typeof providerOrId === 'string' ? this.getById(providerOrId) : providerOrId;
    if (!provider) return;
    const count = this.counts.get(provider);
    if (count === undefined) return;
    if (count > 1) {
      this.counts.set(provider, count - 1);
      return;
    }

    this.counts.delete(provider);
    this.keys.delete(provider.id);
    if (provider.alias !== undefined) this.keys.delete(provider.alias);
    const idx = this._providers.indexOf(provider);
    if (idx >= 0) this._providers.splice(idx, 1);
  }

  has(idOrAlias: string) {
    return this.keys.has(idOrAlias);
  }

  /** Resolves by `id`, then by `alias`. */
  resolve(idOrAlias: string): P | undefined {
    return this.keys.get(idOrAlias);
  }

  private getById(id: string): P | undefined {
    const p = this.keys.get(id);
    return p && p.id === id ? p : undefined;
  }
}
