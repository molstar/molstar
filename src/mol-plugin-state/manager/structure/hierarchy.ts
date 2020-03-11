/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '../../../mol-plugin/context';
import { StructureHierarchy, buildStructureHierarchy, ModelRef, StructureComponentRef, StructureRef, HierarchyRef } from './hierarchy-state';
import { PluginComponent } from '../../component';

interface StructureHierarchyManagerState {
    hierarchy: StructureHierarchy,
    currentModels: ReadonlyArray<ModelRef>,
    currentStructures: ReadonlyArray<StructureRef>
}

export class StructureHierarchyManager extends PluginComponent<StructureHierarchyManagerState> {
    readonly behaviors = {
        current: this.ev.behavior({ hierarchy: this.state.hierarchy, models: this.state.currentModels, structures: this.state.currentStructures })
    }

    private _componentGroups: ReturnType<typeof StructureHierarchyManager['getComponentGroups']> | undefined = void 0;

    get componentGroups() {
        if (this._componentGroups) return this._componentGroups;
        this._componentGroups = StructureHierarchyManager.getComponentGroups(this.state.currentStructures);
        return this._componentGroups;
    }

    private syncCurrentModels(hierarchy: StructureHierarchy): ModelRef[] {
        const current = this.state.currentModels;
        if (current.length === 0) {
            return hierarchy.trajectories[0]?.models || [];
        }

        const newCurrent: ModelRef[] = [];
        for (const c of current) {
            const ref = hierarchy.refs.get(c.cell.transform.ref) as ModelRef;
            if (!ref) continue;
            newCurrent.push(ref);
        }

        if (newCurrent.length === 0) {
            return hierarchy.trajectories[0]?.models || [];
        }

        return newCurrent;
    }

    private syncCurrentStructures(hierarchy: StructureHierarchy, currentModels: ModelRef[]): StructureRef[] {
        const current = this.state.currentStructures;
        if (current.length === 0) {
            return Array.prototype.concat.apply([], currentModels.map(m => m.structures));
        }

        const newCurrent: StructureRef[] = [];
        for (const c of current) {
            const ref = hierarchy.refs.get(c.cell.transform.ref) as StructureRef;
            if (!ref) continue;
            newCurrent.push(ref);
        }

        if (newCurrent.length === 0 && currentModels.length > 0) {
            return Array.prototype.concat.apply([], currentModels.map(m => m.structures));
        }

        return newCurrent;
    }

    private sync() {
        const update = buildStructureHierarchy(this.plugin.state.dataState, this.state.hierarchy);
        if (update.added.length === 0 && update.updated.length === 0 && update.removed.length === 0) {
            return;
        }
        this._componentGroups = void 0;

        const currentModels = this.syncCurrentModels(update.hierarchy);
        const currentStructures = this.syncCurrentStructures(update.hierarchy, currentModels);
        this.updateState({ hierarchy: update.hierarchy, currentModels: currentModels, currentStructures: currentStructures });

        this.behaviors.current.next({ hierarchy: update.hierarchy, models: currentModels, structures: currentStructures });
    }

    remove(refs: HierarchyRef[]) {
        if (refs.length === 0) return;
        const deletes = this.plugin.state.dataState.build();
        for (const r of refs) deletes.delete(r.cell.transform.ref);
        return this.plugin.runTask(this.plugin.state.dataState.updateTree(deletes));
    }

    constructor(private plugin: PluginContext) {
        super({
            hierarchy: StructureHierarchy(),
            currentModels: [],
            currentStructures: []
        });

        plugin.state.dataState.events.changed.subscribe(e => {
            if (e.inTransaction || plugin.behaviors.state.isAnimating.value) return;
            this.sync();
        });

        plugin.behaviors.state.isAnimating.subscribe(isAnimating => {
            if (!isAnimating && !plugin.behaviors.state.isUpdating.value) this.sync();
        });
    }
}

export namespace StructureHierarchyManager {
    export function getComponentGroups(structures: ReadonlyArray<StructureRef>): StructureComponentRef[][] {
        if (!structures.length) return [];
        if (structures.length === 1) return structures[0].components.map(c => [c]);

        const groups: StructureComponentRef[][] = [];
        const map = new Map<string, StructureComponentRef[]>();

        for (const s of structures) {
            for (const c of s.components) {
                const key = c.key;
                if (!key) continue;

                let component = map.get(key);
                if (!component) {
                    component = [];
                    map.set(key, component);
                    groups.push(component);
                }
                component.push(c);
            }
        }

        return groups;
    }
}