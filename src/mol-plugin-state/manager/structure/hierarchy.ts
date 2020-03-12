/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '../../../mol-plugin/context';
import { StructureHierarchy, buildStructureHierarchy, ModelRef, StructureComponentRef, StructureRef, HierarchyRef, TrajectoryRef } from './hierarchy-state';
import { PluginComponent } from '../../component';

interface StructureHierarchyManagerState {
    hierarchy: StructureHierarchy,
    currentTrajectories: ReadonlyArray<TrajectoryRef>,
    currentModels: ReadonlyArray<ModelRef>,
    currentStructures: ReadonlyArray<StructureRef>
}

export class StructureHierarchyManager extends PluginComponent<StructureHierarchyManagerState> {
    readonly behaviors = {
        current: this.ev.behavior({
            hierarchy: this.state.hierarchy,
            trajectories: this.state.currentTrajectories,
            models: this.state.currentModels,
            structures: this.state.currentStructures
        })
    }

    private _currentComponentGroups: ReturnType<typeof StructureHierarchyManager['getComponentGroups']> | undefined = void 0;

    get currentComponentGroups() {
        if (this._currentComponentGroups) return this._currentComponentGroups;
        this._currentComponentGroups = StructureHierarchyManager.getComponentGroups(this.state.currentStructures);
        return this._currentComponentGroups;
    }

    private syncCurrentTrajectories(hierarchy: StructureHierarchy): TrajectoryRef[] {
        const current = this.state.currentTrajectories;
        if (current.length === 0) return hierarchy.trajectories.length > 0 ? [hierarchy.trajectories[0]] : [];

        const newCurrent: TrajectoryRef[] = [];
        for (const c of current) {
            const ref = hierarchy.refs.get(c.cell.transform.ref) as TrajectoryRef;
            if (ref) newCurrent.push(ref);
        }

        if (newCurrent.length === 0) return hierarchy.trajectories.length > 0 ? [hierarchy.trajectories[0]] : [];
        return newCurrent;
    }

    private syncCurrentModels(hierarchy: StructureHierarchy, currentTrajectories: TrajectoryRef[]): ModelRef[] {
        const current = this.state.currentModels;
        if (current.length === 0) return currentTrajectories[0]?.models || [];

        const newCurrent: ModelRef[] = [];
        for (const c of current) {
            const ref = hierarchy.refs.get(c.cell.transform.ref) as ModelRef;
            if (ref) newCurrent.push(ref);
        }

        if (newCurrent.length === 0) return currentTrajectories[0]?.models || [];
        return newCurrent;
    }

    private syncCurrentStructures(hierarchy: StructureHierarchy, currentModels: ModelRef[]): StructureRef[] {
        const current = this.state.currentStructures;
        if (current.length === 0) return Array.prototype.concat.apply([], currentModels.map(m => m.structures));

        const newCurrent: StructureRef[] = [];
        for (const c of current) {
            const ref = hierarchy.refs.get(c.cell.transform.ref) as StructureRef;
            if (ref) newCurrent.push(ref);
        }

        if (newCurrent.length === 0 && currentModels.length > 0) return Array.prototype.concat.apply([], currentModels.map(m => m.structures));
        return newCurrent;
    }

    private sync() {
        const update = buildStructureHierarchy(this.plugin.state.dataState, this.state.hierarchy);
        if (update.added.length === 0 && update.updated.length === 0 && update.removed.length === 0) {
            return;
        }
        this._currentComponentGroups = void 0;

        const currentTrajectories = this.syncCurrentTrajectories(update.hierarchy);
        const currentModels = this.syncCurrentModels(update.hierarchy, currentTrajectories);
        const currentStructures = this.syncCurrentStructures(update.hierarchy, currentModels);
        console.log(currentTrajectories, currentModels, currentStructures);
        this.updateState({ hierarchy: update.hierarchy, currentModels, currentStructures, currentTrajectories });

        this.behaviors.current.next({
            hierarchy: update.hierarchy,
            trajectories: currentTrajectories,
            models: currentModels,
            structures: currentStructures
        });
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
            currentTrajectories: [],
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