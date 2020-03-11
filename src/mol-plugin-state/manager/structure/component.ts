/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureRef } from './hierarchy-state'
import { StructureRepresentationProvider } from '../../builder/structure/provider';
import { PluginContext } from '../../../mol-plugin/context';

export { StructureComponentManager }

class StructureComponentManager {
    applyPreset<P = any, S = {}>(structures: StructureRef[], provider: StructureRepresentationProvider<P, S>, params?: P): Promise<any>  {
        return this.plugin.runTask(this.dateState.transaction(async () => {
            await this.removeComponents(structures);
            for (const s of structures) {
                await this.plugin.builders.structure.representation.structurePreset(s.cell, provider, params);
            }
        }));
    }

    clear(structures: StructureRef[]) {
        return this.removeComponents(structures);
    }

    modify(structures: StructureRef[], action: StructureComponentManager.ModifyAction) {

    }

    private get dateState() {
        return this.plugin.state.dataState;
    }

    private removeComponents(structures: StructureRef[]) {
        const deletes = this.dateState.build();
        for (const s of structures) {
            for (const c of s.components) {
                deletes.delete(c.cell.transform.ref);
            }
            if (s.currentFocus) {
                if (s.currentFocus.focus) deletes.delete(s.currentFocus.focus.cell.transform.ref);
                if (s.currentFocus.surroundings) deletes.delete(s.currentFocus.surroundings.cell.transform.ref);
            }
        }
        return this.plugin.runTask(this.dateState.updateTree(deletes));
    }

    constructor(public plugin: PluginContext) {

    }
}

namespace StructureComponentManager {

    export function getModifyParams() {
        return 0 as any;
    }

    export type ModifyAction = 
        | { kind: 'add', label: string, representationType?: string }
        | { kind: 'merge', type: { kind: 'intersecting', key: string } | { kind: 'component', key: string } }
        | { kind: 'subtract', type: { kind: 'all' } | { kind: 'component', key: string } }
        | { kind: 'color', representationType?: string }
}
