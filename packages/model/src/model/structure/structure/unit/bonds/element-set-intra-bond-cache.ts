/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StructureElement } from '../../element.js';
import type { IntraUnitBonds } from './data.js';
import { SortedArray } from '@molstar/core/data/int';
import type { Model } from '../../../model.js';
import { LRUCache } from '@molstar/core/util/lru-cache';

type Entry = [StructureElement.Set, IntraUnitBonds];

export class ElementSetIntraBondCache {
  private data = new Map<number, Entry[]>();
  /** The entries in order of use, least recent first: past `MaxEntries` the least recent goes (perf-audit-2 05/P11) */
  private order = new Set<Entry>();

  get(xs: StructureElement.Set): IntraUnitBonds | undefined {
    const hash = SortedArray.hashCode(xs);
    if (!this.data.has(hash)) return void 0;
    for (const e of this.data.get(hash)!) {
      if (SortedArray.areEqual(xs, e[0])) {
        this.order.delete(e);
        this.order.add(e);
        return e[1];
      }
    }
  }

  set(xs: StructureElement.Set, bonds: IntraUnitBonds) {
    const hash = SortedArray.hashCode(xs);
    let es = this.data.get(hash);
    if (es) {
      for (const e of es) {
        if (SortedArray.areEqual(xs, e[0])) {
          e[1] = bonds;
          this.order.delete(e);
          this.order.add(e);
          return;
        }
      }
    } else {
      this.data.set(hash, (es = []));
    }
    const e: Entry = [xs, bonds];
    es.push(e);
    this.order.add(e);
    if (this.order.size > ElementSetIntraBondCache.MaxEntries) this.evict();
  }

  private evict() {
    const e = this.order.values().next().value!;
    this.order.delete(e);
    const hash = SortedArray.hashCode(e[0]);
    const es = this.data.get(hash)!;
    es.splice(es.indexOf(e), 1);
    if (es.length === 0) this.data.delete(hash);
  }

  /**
   * How many element sets a cache keeps bonds for: each selection or component of an index-pair model adds one, and
   * units hold the bonds they have, so an evicted set costs only a recompute.
   */
  static readonly MaxEntries = 64;

  /**
   * The model's cache for the bonds of one provider instance's key (`BondProvider.Instance.cacheKey`): two structures
   * of a model with different providers never share bonds. The default key is the model's original cache, kept for
   * the model's life; of the others only the `MaxProviderKeys` used last are kept, as a key per param value (a cutoff
   * slider) would otherwise keep every value's bonds. Units hold the bonds they have, so a dropped key costs only a
   * recompute.
   */
  static get(model: Model, key = 'default'): ElementSetIntraBondCache {
    const data = model._dynamicPropertyData;
    if (key === 'default') {
      if (!data.ElementSetIntraBondCache) data.ElementSetIntraBondCache = new ElementSetIntraBondCache();
      return data.ElementSetIntraBondCache;
    }
    let byKey: LRUCache<ElementSetIntraBondCache> | undefined = data.ElementSetIntraBondCacheByKey;
    if (!byKey) data.ElementSetIntraBondCacheByKey = byKey = LRUCache.create(ElementSetIntraBondCache.MaxProviderKeys);
    let cache = LRUCache.get(byKey, key);
    if (!cache) LRUCache.set(byKey, key, (cache = new ElementSetIntraBondCache()));
    return cache;
  }

  /** How many provider keys a model keeps caches for besides the default. */
  static readonly MaxProviderKeys = 4;
}
