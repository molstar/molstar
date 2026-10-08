/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { FileNameInfo } from '@molstar/core/util/file-info';
import type { PluginStateObject } from '../objects.js';
import { DataFormatProvider } from './provider.js';

const warnedRenames = new WeakMap<object, Set<string>>();

function warnLegacyRename(provider: DataFormatProvider.Unnamed, name: string) {
  let names = warnedRenames.get(provider);
  if (!names) {
    names = new Set();
    warnedRenames.set(provider, names);
  }
  if (names.has(name)) return;
  names.add(name);
  console.warn(
    `DataFormatRegistry.add(name, provider) with a name that differs from provider.name ('${provider.name ?? ''}') is deprecated. Register DataFormatProvider.withName(provider, '${name}') instead.`,
  );
}

export class DataFormatRegistry {
  private _list: { name: string; provider: DataFormatProvider }[] = [];
  private _map = new Map<string, { provider: DataFormatProvider; count: number }>();
  private _extensions: Set<string> | undefined = undefined;
  private _binaryExtensions: Set<string> | undefined = undefined;
  private _options: [name: string, label: string, category: string][] | undefined = undefined;
  private _autoOrder: { name: string; provider: DataFormatProvider }[] | undefined = undefined;

  get types(): [name: string, label: string][] {
    return this._list.map((e) => [e.name, e.provider.label] as [name: string, label: string]);
  }

  get extensions() {
    if (this._extensions) return this._extensions;
    const extensions = new Set<string>();
    this._list.forEach(({ provider }) => {
      provider.stringExtensions?.forEach((ext) => extensions.add(ext));
      provider.binaryExtensions?.forEach((ext) => extensions.add(ext));
    });
    this._extensions = extensions;
    return extensions;
  }

  get binaryExtensions() {
    if (this._binaryExtensions) return this._binaryExtensions;
    const binaryExtensions = new Set<string>();
    this._list.forEach(({ provider }) => provider.binaryExtensions?.forEach((ext) => binaryExtensions.add(ext)));
    this._binaryExtensions = binaryExtensions;
    return binaryExtensions;
  }

  get options() {
    if (this._options) return this._options;
    const options: [name: string, label: string, category: string][] = [];
    this._list.forEach(({ name, provider }) => options.push([name, provider.label, provider.category || '']));
    this._options = options;
    return options;
  }

  /**
   * Providers ordered for `auto()`: higher `priority` first, ties broken by registration
   * order (stable sort). Explicit priority makes auto-detection independent of the order
   * providers happen to be registered in.
   */
  get autoOrder() {
    if (this._autoOrder) return this._autoOrder;
    this._autoOrder = this._list
      .map((entry, index) => ({ entry, index }))
      .sort((a, b) => (b.entry.provider.priority ?? 0) - (a.entry.provider.priority ?? 0) || a.index - b.index)
      .map(({ entry }) => entry);
    return this._autoOrder;
  }

  private _invalidate() {
    this._extensions = undefined;
    this._binaryExtensions = undefined;
    this._options = undefined;
    this._autoOrder = undefined;
  }

  /**
   * Returns the message `add` would throw for `provider`, or `undefined` when it can be added.
   * Does not change the registry.
   */
  findConflict(provider: DataFormatProvider): string | undefined {
    const existing = this._map.get(provider.name);
    if (existing && existing.provider !== provider) {
      return `Data format registry: '${provider.name}' is already registered by a different provider.`;
    }
    return undefined;
  }

  /**
   * Registers a provider under `provider.name`. Adding the same object again increments its count,
   * a different object under an existing name throws.
   *
   * `add(name, provider)` is deprecated: when `name` differs from `provider.name` it registers
   * `DataFormatProvider.withName(provider, name)` instead.
   */
  add(provider: DataFormatProvider): void;
  /** @deprecated Use `add(provider)`, and `DataFormatProvider.withName` to register a provider under another name. */
  add(name: string, provider: DataFormatProvider.Unnamed): void;
  add(nameOrProvider: string | DataFormatProvider, maybeProvider?: DataFormatProvider.Unnamed) {
    let provider: DataFormatProvider;
    if (typeof nameOrProvider === 'string') {
      const name = nameOrProvider;
      const source = maybeProvider!;
      provider = DataFormatProvider.withName(source, name);
      if (provider !== source) warnLegacyRename(source, name);
    } else {
      provider = nameOrProvider;
    }

    const conflict = this.findConflict(provider);
    if (conflict) throw new Error(conflict);

    const existing = this._map.get(provider.name);
    if (existing) {
      existing.count++;
      return;
    }
    this._invalidate();
    this._map.set(provider.name, { provider, count: 1 });
    this._list.push({ name: provider.name, provider });
  }

  /** Decrements the count of a provider (given by identity or name) and removes it at zero. Unknown is a no-op. */
  remove(providerOrName: DataFormatProvider | string) {
    const name = typeof providerOrName === 'string' ? providerOrName : providerOrName.name;
    const existing = this._map.get(name);
    if (!existing) return;
    if (typeof providerOrName !== 'string' && existing.provider !== providerOrName) return;
    if (--existing.count > 0) return;

    this._invalidate();
    this._map.delete(name);
    const index = this._list.findIndex((e) => e.name === name);
    if (index >= 0) this._list.splice(index, 1);
  }

  /** Drops every provider and count. */
  clear() {
    this._invalidate();
    this._map.clear();
    this._list = [];
  }

  has(name: string) {
    return this._map.has(name);
  }

  auto(info: FileNameInfo, dataStateObject: PluginStateObject.Data.Binary | PluginStateObject.Data.String) {
    const list = this.autoOrder;
    for (let i = 0, il = list.length; i < il; ++i) {
      const p = list[i].provider;

      let hasExt = false;
      if (p.binaryExtensions?.includes(info.ext)) hasExt = true;
      else if (p.stringExtensions?.includes(info.ext)) hasExt = true;
      if (hasExt && (!p.isApplicable || p.isApplicable(info, dataStateObject.data))) return p;
    }
    return;
  }

  get(name: string): DataFormatProvider | undefined {
    const entry = this._map.get(name);
    if (!entry) throw new Error(unregisteredFormatMessage(name, this._list.length === 0));
    return entry.provider;
  }

  get list() {
    return this._list;
  }
}

/** The most common cause of an empty format registry is a 5.x-style spec without `registry`, so say how to fix it. */
export function unregisteredFormatMessage(name: string, registryIsEmpty: boolean) {
  const message = `Data format '${name}' is not registered in this plugin.`;
  if (!registryIsEmpty) return message;
  return `${message} No data formats are registered: plugin registries start empty, so add format entries to \`spec.registry\` (for the full built-in set, use \`DefaultRegistry\` from '@molstar/plugin/default-registry' or \`DefaultPluginSpec()\`).`;
}
