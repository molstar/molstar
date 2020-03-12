/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '../../../mol-plugin/context';
import { StructureHierarchy, buildStructureHierarchy, ModelRef, StructureComponentRef, StructureRef, HierarchyRef, TrajectoryRef } from './hierarchy-state';
import { PluginComponent } from '../../component';
import { StateTransform } from '../../../mol-state';
import { SetUtils } from '../../../mol-util/set';

interface StructureHierarchyManagerState {
    hierarchy: StructureHierarchy,
    current: {
        trajectories: ReadonlyArray<TrajectoryRef>,
        models: ReadonlyArray<ModelRef>,
        structures: ReadonlyArray<StructureRef>
    }
}

export class StructureHierarchyManager extends PluginComponent<StructureHierarchyManagerState> {
    readonly behaviors = {
        current: this.ev.behavior({
            hierarchy: this.state.hierarchy,
            trajectories: this.state.current.trajectories,
            models: this.state.current.models,
            structures: this.state.current.structures
        })
    }

    private get dataState() {
        return this.plugin.state.dataState;
    }

    private _currentComponentGroups: ReturnType<typeof StructureHierarchyManager['getComponentGroups']> | undefined = void 0;

    get currentComponentGroups() {
        if (this._currentComponentGroups) return this._currentComponentGroups;
        this._currentComponentGroups = StructureHierarchyManager.getComponentGroups(this.state.current.structures);
        return this._currentComponentGroups;
    }

    private _currentSelectionSet: Set<string> | undefined = void 0;
    get currentSeletionSet() {
        if (this._currentSelectionSet) return this._currentSelectionSet;
        this._currentSelectionSet = new Set();
        for (const r of this.state.current.trajectories) this._currentSelectionSet.add(r.cell.transform.ref);
        for (const r of this.state.current.models) this._currentSelectionSet.add(r.cell.transform.ref);
        for (const r of this.state.current.structures) this._currentSelectionSet.add(r.cell.transform.ref);
        return this._currentSelectionSet;
    }

    private syncCurrentTrajectories(hierarchy: StructureHierarchy, map: Map<StateTransform.Ref, HierarchyRef>): TrajectoryRef[] {
        const current = this.state.current.trajectories;
        if (current.length === 0) return hierarchy.trajectories.length > 0 ? [hierarchy.trajectories[0]] : [];

        const newCurrent: TrajectoryRef[] = [];
        for (const c of current) {
            const ref = map.get(c.cell.transform.ref) as TrajectoryRef;
            if (ref) newCurrent.push(ref);
        }

        if (newCurrent.length === 0) return hierarchy.trajectories.length > 0 ? [hierarchy.trajectories[0]] : [];
        return newCurrent;
    }

    private syncCurrentModels(hierarchy: StructureHierarchy, map: Map<StateTransform.Ref, HierarchyRef>, currentTrajectories: TrajectoryRef[]): ModelRef[] {
        const current = this.state.current.models;
        if (current.length === 0) return currentTrajectories[0]?.models || [];

        const newCurrent: ModelRef[] = [];
        for (const c of current) {
            const ref = map.get(c.cell.transform.ref) as ModelRef;
            if (ref) newCurrent.push(ref);
        }

        if (newCurrent.length === 0) return currentTrajectories[0]?.models || [];
        return newCurrent;
    }

    private syncCurrentStructures(map: Map<StateTransform.Ref, HierarchyRef>, currentModels: ModelRef[]): StructureRef[] {
        const current = this.state.current.structures;
        if (current.length === 0) return Array.prototype.concat.apply([], currentModels.map(m => m.structures));

        const newCurrent: StructureRef[] = [];
        for (const c of current) {
            const ref = map.get(c.cell.transform.ref) as StructureRef;
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
        this._currentSelectionSet = void 0;

        const trajectories = this.syncCurrentTrajectories(update.hierarchy, update.hierarchy.refs);
        const models = this.syncCurrentModels(update.hierarchy, update.hierarchy.refs, trajectories);
        const structures = this.syncCurrentStructures(update.hierarchy.refs, models);

        this.updateState({ hierarchy: update.hierarchy, current: { trajectories, models, structures }});
        this.behaviors.current.next({ hierarchy: update.hierarchy, trajectories, models, structures });
    }

    updateCurrent(refs: HierarchyRef[], action: 'add' | 'remove') {
        const hierarchy = this.state.hierarchy;
        const set = action === 'add'
            ? SetUtils.union(this.currentSeletionSet, new Set(refs.map(r => r.cell.transform.ref)))
            : SetUtils.difference(this.currentSeletionSet, new Set(refs.map(r => r.cell.transform.ref)));

        const trajectories = [];
        const models = [];
        const structures = [];

        for (const t of hierarchy.trajectories) {
            if (set.has(t.cell.transform.ref)) trajectories.push(t);
            for (const m of t.models) {
                if (set.has(m.cell.transform.ref)) models.push(m);
                for (const s of m.structures) {
                    if (set.has(s.cell.transform.ref)) structures.push(s);
                }
            }
        }

        this._currentComponentGroups = void 0;
        this._currentSelectionSet = void 0;
     
        this.updateState({ current: { trajectories, models, structures }});
        this.behaviors.current.next({ hierarchy, trajectories, models, structures });
    }
    
    remove(refs: HierarchyRef[]) {
        if (refs.length === 0) return;
        const deletes = this.plugin.state.dataState.build();
        for (const r of refs) deletes.delete(r.cell.transform.ref);
        return this.plugin.runTask(this.plugin.state.dataState.updateTree(deletes));
    }

    createAllModels(trajectory: TrajectoryRef) {
        return this.plugin.dataTransaction(async () => {
            if (trajectory.models.length > 0) {
                await this.clearTrajectory(trajectory);
            }

            const tr = trajectory.cell.obj?.data!;
            for (let i = 0; i < tr.length; i++) {
                const model = await this.plugin.builders.structure.createModel(trajectory.cell, { modelIndex: i });
                const structure = await this.plugin.builders.structure.createStructure(model, { name: 'deposited', params: { } });
                await this.plugin.builders.structure.representation.structurePreset(structure, 'auto');
            }
        })
    }

    private clearTrajectory(trajectory: TrajectoryRef) {
        const builder = this.dataState.build();
        for (const m of trajectory.models) {
            builder.delete(m.cell);
        }
        return this.plugin.runTask(this.dataState.updateTree(builder));
    }

    constructor(private plugin: PluginContext) {
        super({
            hierarchy: StructureHierarchy(),
            current: { trajectories: [], models: [], structures: [] }
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