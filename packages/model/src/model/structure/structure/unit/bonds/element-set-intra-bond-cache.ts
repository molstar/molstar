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

export class ElementSetIntraBondCache {
  private data = new Map<number, [StructureElement.Set, IntraUnitBonds][]>();

  get(xs: StructureElement.Set): IntraUnitBonds | undefined {
    const hash = SortedArray.hashCode(xs);
    if (!this.data.has(hash)) return void 0;
    for (const [s, b] of this.data.get(hash)!) {
      if (SortedArray.areEqual(xs, s)) return b;
    }
  }

  set(xs: StructureElement.Set, bonds: IntraUnitBonds) {
    const hash = SortedArray.hashCode(xs);
    if (this.data.has(hash)) {
      const es = this.data.get(hash)!;
      for (const e of es) {
        if (SortedArray.areEqual(xs, e[0])) {
          e[1] = bonds;
          return;
        }
      }
      es.push([xs, bonds]);
    } else {
      this.data.set(hash, [[xs, bonds]]);
    }
  }

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
