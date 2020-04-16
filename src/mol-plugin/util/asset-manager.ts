/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { UUID } from '../../mol-util';
import { PluginContext } from '../context';
import { iterableToArray } from '../../mol-data/util';

export { AssetManager };

class AssetManager {
    private files = new Map<UUID, File>()
    private ids = new Map<File, UUID>()

    get list() {
        return iterableToArray(this.ids.entries());
    }

    set(id: UUID, file: File) {
        this.files.set(id, file);
        this.ids.set(file, id);
    }

    remove(id: UUID) {
        if (this.files.has(id)) {
            const file = this.files.get(id)!;
            this.files.delete(id);
            this.ids.delete(file);
        }
    }

    has(id: UUID) {
        return this.files.has(id);
    }

    get(id: UUID) {
        return this.files.get(id);
    }

    /** For use with `JSON.stringify` */
    replacer = (key: string, value: any) => {
        if (value instanceof File) {
            const id = this.ids.get(value);
            if (!id) {
                // TODO throw?
                console.warn(`No asset found for '${value.name}'`);
            }
            return id ? AssetManager.Item(id, value.name) : {};
        } else {
            return value;
        }
    }

    constructor(public ctx: PluginContext) {
        ctx.state.data.events.object.removed.subscribe(e => {
            const id = e.obj?.id;
            if (id) ctx.managers.asset.remove(id);
        });
    }
}

namespace AssetManager {
    export type Item = { kind: 'asset-item', id: UUID, name: string };
    export function Item(id: UUID, name: string): Item {
        return { kind: 'asset-item', id, name };
    }

    export function isItem(x?: any): x is Item {
        return !!x && x?.kind === 'asset-item';
    }
}