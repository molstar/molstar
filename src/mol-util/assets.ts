/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import UUID from './uuid';
import { iterableToArray } from '../mol-data/util';
import { ajaxGet, DataType, DataResponse, readFromFile } from './data-source';
import { Task } from '../mol-task';

export { AssetManager, Asset };

type _File = File;
type Asset = Asset.Url | Asset.File

namespace Asset {
    export type Url = { url: string, title?: string, body?: string }
    export type File = { id: UUID, name: string, file?: _File }

    export function Url(url: string, options?: { body?: string, title?: string }): Url {
        return { url, ...options };
    }

    export function File(file: _File): File {
        return { id: UUID.create22(), name: file.name, file };
    }

    export function isUrl(x: Asset): x is Url {
        return !!x && !!(x as any).url;
    }

    export function isFile(x: Asset): x is File {
        return !!x && !!(x as any).id;
    }
}

class AssetManager {
    private _assets = new Map<string, { asset: Asset, file: File }>();

    get assets() {
        return iterableToArray(this._assets.values());
    }

    set(asset: Asset, file: File) {
        if (Asset.isUrl(asset)) {
            this._assets.set(getUrlKey(asset), { asset, file });
        } else {
            this._assets.set(asset.id, { asset, file });
        }
    }

    resolve<T extends DataType>(asset: Asset, type: T, store = true): Task<DataResponse<T>> {
        if (Asset.isUrl(asset)) {
            const key = getUrlKey(asset);
            if (this._assets.has(key)) {
                return readFromFile(this._assets.get(key)!.file, type);
            }

            if (!store) {
                return ajaxGet({ ...asset, type });
            }

            return Task.create(`Download ${asset.title || asset.url}`, async ctx => {
                const data = await ajaxGet({ ...asset, type: 'binary' }).runInContext(ctx);
                const file = new File([data], 'raw-data');
                this._assets.set(key, { asset, file });
                return await readFromFile(file, type).runInContext(ctx);
            });
        } else {
            if (this._assets.has(asset.id)) return readFromFile(this._assets.get(asset.id)!.file, type);
            if (!(asset.file instanceof File)) {
                return Task.fail('Resolve asset', `Cannot resolve file asset '${asset.name}' (${asset.id})`);
            }
            if (store) {
                this._assets.set(asset.id, { asset, file: asset.file });
            }
            return readFromFile(asset.file, type);
        }
    }

    release(asset: Asset) {
        if (Asset.isFile(asset)) {
            this._assets.delete(asset.id);
        } else {
            this._assets.delete(getUrlKey(asset));
        }
    }
}

function getUrlKey(asset: Asset.Url) {
    return asset.body ? `${asset.url}_${asset.body || ''}` : asset.url;
}