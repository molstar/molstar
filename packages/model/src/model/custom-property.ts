/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { CifWriter } from '@molstar/io/writer/cif';
import type { CifExportContext } from './structure/export/mmcif.js';
import type { QuerySymbolRuntime } from '@molstar/model/script/runtime/query/compiler';
import { UUID } from '@molstar/core/util';

export { CustomPropertyDescriptor, CustomProperties };

interface CustomPropertyDescriptor<
  ExportCtx = CifExportContext,
  Symbols extends { [name: string]: QuerySymbolRuntime } = {},
> {
  readonly name: string;

  cifExport?: {
    // Prefix enforced during export.
    prefix: string;
    context?: (ctx: CifExportContext) => ExportCtx | undefined;
    categories: CifWriter.Category<ExportCtx>[];
  };

  // TODO: add aliases when lisp-like mol-script is done
  symbols?: Symbols;
}

function CustomPropertyDescriptor<Ctx, Desc extends CustomPropertyDescriptor<Ctx>>(desc: Desc) {
  return desc;
}

namespace CustomPropertyDescriptor {
  export function getUUID(prop: CustomPropertyDescriptor): UUID {
    if (!(prop as any).__key) {
      (prop as any).__key = UUID.create22();
    }
    return (prop as any).__key;
  }
}

/**
 * Anything with a dispose method, used to despose of data assets or webgl resources
 */
type Asset = { dispose: () => void };

class CustomProperties {
  private _list: CustomPropertyDescriptor[] = [];
  private _set = new Set<CustomPropertyDescriptor>();
  private _refs = new Map<CustomPropertyDescriptor, number>();
  /** desc → owner → assets; the owner is `undefined` unless the objects sharing these properties own them apart */
  private _assets = new Map<CustomPropertyDescriptor, Map<string | undefined, Asset[]>>();

  get all(): ReadonlyArray<CustomPropertyDescriptor> {
    return this._list;
  }

  add(desc: CustomPropertyDescriptor<any>) {
    if (this._set.has(desc)) return;

    this._list.push(desc);
    this._set.add(desc);
  }

  reference(desc: CustomPropertyDescriptor<any>, add: boolean) {
    let refs = this._refs.get(desc) || 0;
    refs += add ? 1 : -1;
    this._refs.set(desc, Math.max(refs, 0));
  }

  hasReference(desc: CustomPropertyDescriptor<any>) {
    return (this._refs.get(desc) || 0) > 0;
  }

  has(desc: CustomPropertyDescriptor<any>): boolean {
    return this._set.has(desc);
  }

  /**
   * Sets assets for a prop, disposes of existing assets for that prop. `owner` keeps the assets of one of the objects
   * sharing these properties apart, e.g. the model id for dynamic data: models with the same topology share
   * `customProperties` but not their dynamic property data.
   */
  assets(desc: CustomPropertyDescriptor<any>, assets?: Asset[], owner?: string) {
    let byOwner = this._assets.get(desc);
    const prevAssets = byOwner?.get(owner);
    if (prevAssets) {
      for (const a of prevAssets) a.dispose();
    }
    if (assets) {
      if (!byOwner) {
        byOwner = new Map();
        this._assets.set(desc, byOwner);
      }
      byOwner.set(owner, assets);
    } else if (byOwner) {
      byOwner.delete(owner);
      if (byOwner.size === 0) this._assets.delete(desc);
    }
  }

  /** Disposes of all assets of all props (and forgets them, so a later `assets()` does not dispose them again) */
  dispose() {
    this._assets.forEach((byOwner) => {
      byOwner.forEach((assets) => {
        for (const a of assets) a.dispose();
      });
    });
    this._assets.clear();
  }
}
