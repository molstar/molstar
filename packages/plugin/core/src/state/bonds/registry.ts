/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Model } from '@molstar/model/model/structure';
import { BondProvider } from '@molstar/model/model/structure/structure/unit/bonds/provider';

/** A stored choice of bond provider: what the state tree holds (bond-providers plan §4). */
export interface BondProviderChoice {
  name: string;
  params: any;
}

/**
 * The bond providers a plugin offers, keyed by name and counted like the other registries (plugin-composition §4.2):
 * the same object again counts up, a different object under a taken name throws, removal counts down. Providers may
 * be added and removed at any time (bond-providers plan §3); existing structures keep the instances they hold.
 *
 * `'default'` (`BondProvider.Default`, the built-in bonds) is always there and is not an entry: it cannot be added or
 * removed, it is the first option, and every name that does not resolve falls back to it.
 */
export class BondProviderRegistry {
  private _list: BondProvider[] = [];
  private _map = new Map<string, { provider: BondProvider; count: number }>();

  /** The registered providers, in registration order, without `'default'`. */
  get list(): readonly BondProvider[] {
    return this._list;
  }

  /** Returns the message `add` would throw for `provider`, or `undefined` when it can be added. */
  findConflict(provider: BondProvider): string | undefined {
    if (provider.name === BondProvider.DefaultName) {
      return `Bond provider registry: '${BondProvider.DefaultName}' is the built-in provider and cannot be registered.`;
    }
    const existing = this._map.get(provider.name);
    if (existing && existing.provider !== provider) {
      return `Bond provider registry: '${provider.name}' is already registered by a different provider.`;
    }
    return undefined;
  }

  add(provider: BondProvider) {
    const conflict = this.findConflict(provider);
    if (conflict) throw new Error(conflict);
    const existing = this._map.get(provider.name);
    if (existing) {
      existing.count++;
      return;
    }
    this._map.set(provider.name, { provider, count: 1 });
    this._list.push(provider);
  }

  /** Decrements the count of a provider (given by identity or name) and removes it at zero. Unknown is a no-op. */
  remove(providerOrName: BondProvider | string) {
    const name = typeof providerOrName === 'string' ? providerOrName : providerOrName.name;
    const existing = this._map.get(name);
    if (!existing) return;
    if (typeof providerOrName !== 'string' && existing.provider !== providerOrName) return;
    if (--existing.count > 0) return;
    this._map.delete(name);
    this._list.splice(this._list.indexOf(existing.provider), 1);
  }

  clear() {
    this._map.clear();
    this._list = [];
  }

  /** `'default'` is always there. */
  has(name: string) {
    return name === BondProvider.DefaultName || this._map.has(name);
  }

  /** The provider of `name`: `BondProvider.Default` for `'default'`, `undefined` when not registered. */
  get(name: string): BondProvider | undefined {
    if (name === BondProvider.DefaultName) return BondProvider.Default;
    return this._map.get(name)?.provider;
  }

  /** `'default'`, then the registered providers applicable to `model` (all of them without a model). */
  getApplicable(model?: Model): BondProvider[] {
    const ret: BondProvider[] = [BondProvider.Default];
    for (const p of this._list) {
      if (!model || !p.isApplicable || p.isApplicable(model)) ret.push(p);
    }
    return ret;
  }

  options(model?: Model): [name: string, label: string][] {
    return this.getApplicable(model).map((p) => [p.name, p.label]);
  }

  /**
   * The param of a choice of provider: a `PD.Mapped` over the providers applicable to `model`, `'default'` first, each
   * with its own params for the model. Normalization against it maps an unregistered or non-applicable name to
   * `'default'` (bond-providers plan §4). Never throws, with or without a model.
   */
  getParam(model?: Model, info?: PD.Info) {
    return PD.Mapped<any>(
      BondProvider.DefaultName,
      this.options(model),
      (name) => PD.Group(this.get(name)?.getParams(model) ?? {}),
      { label: 'Bonds', description: 'Where the bonds come from', ...info },
    );
  }

  /**
   * The instance for a stored choice: the default instance for no choice, `'default'`, an unregistered name, or a
   * provider not applicable to `model`; else the provider's instance for its params, normalized against its params for
   * the model.
   */
  resolve(choice: BondProviderChoice | undefined, model?: Model): BondProvider.Instance {
    if (!choice || choice.name === BondProvider.DefaultName) return BondProvider.DefaultInstance;
    const provider = this._map.get(choice.name)?.provider;
    if (!provider) return BondProvider.DefaultInstance;
    if (model && provider.isApplicable && !provider.isApplicable(model)) return BondProvider.DefaultInstance;
    const params = PD.normalizeParams(provider.getParams(model), choice.params ?? {}, 'all');
    return provider.create(params);
  }
}

/** Says when a stored choice names a provider this plugin does not have; the built-in bonds are used instead. */
export function unregisteredBondProviderMessage(name: string) {
  return `Bond provider '${name}' is not registered in this plugin; the default bonds are used.`;
}
