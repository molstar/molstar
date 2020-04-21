/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { setSubtreeVisibility } from '../../../mol-plugin/behavior/static/state';
import { PluginCommands } from '../../../mol-plugin/commands';
import { PluginContext } from '../../../mol-plugin/context';
import { StateTransform, StateTree } from '../../../mol-state';
import { SetUtils } from '../../../mol-util/set';
import { TrajectoryHierarchyPresetProvider } from '../../builder/structure/hierarchy-preset';
import { PluginComponent } from '../../component';
import { buildStructureHierarchy, StructureHierarchyRef, ModelRef, StructureComponentRef, StructureHierarchy, StructureRef, TrajectoryRef } from './hierarchy-state';

export class StructureHierarchyManager extends PluginComponent {
    private state = {
        syncedTree: this.dataState.tree,
        notified: false,

        hierarchy: StructureHierarchy(),
        selection: {
            trajectories: [] as ReadonlyArray<TrajectoryRef>,
            models: [] as ReadonlyArray<ModelRef>,
            structures: [] as ReadonlyArray<StructureRef>
        }
    }

    readonly behaviors = {
        selection: this.ev.behavior({
            hierarchy: this.current,
            trajectories: this.selection.trajectories,
            models: this.selection.models,
            structures: this.selection.structures
        })
    }

    private get dataState() {
        return this.plugin.state.data;
    }

    private _currentComponentGroups: ReturnType<typeof StructureHierarchyManager['getComponentGroups']> | undefined = void 0;

    get currentComponentGroups() {
        if (this._currentComponentGroups) return this._currentComponentGroups;
        this._currentComponentGroups = StructureHierarchyManager.getComponentGroups(this.selection.structures);
        return this._currentComponentGroups;
    }

    private _currentSelectionSet: Set<StateTransform.Ref> | undefined = void 0;
    get seletionSet() {
        if (this._currentSelectionSet) return this._currentSelectionSet;
        this._currentSelectionSet = new Set();
        for (const r of this.selection.trajectories) this._currentSelectionSet.add(r.cell.transform.ref);
        for (const r of this.selection.models) this._currentSelectionSet.add(r.cell.transform.ref);
        for (const r of this.selection.structures) this._currentSelectionSet.add(r.cell.transform.ref);
        return this._currentSelectionSet;
    }

    get current() {
        this.sync(false);
        return this.state.hierarchy;
    }

    get selection() {
        this.sync(false);
        return this.state.selection;
    }

    getStructuresWithSelection() {
        const xs = this.plugin.managers.structure.hierarchy.current.structures;
        const ret: StructureRef[] = [];
        for (const s of xs) {
            if (this.plugin.managers.structure.selection.structureHasSelection(s)) {
                ret.push(s);
            }
        }
        return ret;
    }

    private syncCurrent<T extends StructureHierarchyRef>(all: ReadonlyArray<T>, added: Set<StateTransform.Ref>): T[] {
        const current = this.seletionSet;
        const newCurrent: T[] = [];

        for (const r of all) {
            const ref = r.cell.transform.ref;
            if (current.has(ref) || added.has(ref)) newCurrent.push(r);
        }

        if (newCurrent.length === 0) return all.length > 0 ? [all[0]] : [];
        return newCurrent;
    }

    private sync(notify: boolean) {
        if (!notify && this.dataState.inUpdate) return;

        if (this.state.syncedTree === this.dataState.tree) {
            if (notify && !this.state.notified) {
                this.state.notified = true;
                this.behaviors.selection.next({ hierarchy: this.state.hierarchy, ...this.state.selection });
            }

            return;
        }

        this.state.syncedTree = this.dataState.tree;

        const update = buildStructureHierarchy(this.plugin.state.data, this.current);
        if (!update.changed) {
            return;
        }

        const { hierarchy } = update;
        const trajectories = this.syncCurrent(hierarchy.trajectories, update.added);
        const models = this.syncCurrent(hierarchy.models, update.added);
        const structures = this.syncCurrent(hierarchy.structures, update.added);

        this._currentComponentGroups = void 0;
        this._currentSelectionSet = void 0;

        this.state.hierarchy = hierarchy;
        this.state.selection.trajectories = trajectories;
        this.state.selection.models = models;
        this.state.selection.structures = structures;

        if (notify) {
            this.state.notified = true;
            this.behaviors.selection.next({ hierarchy, trajectories, models, structures });
        } else {
            this.state.notified = false;
        }
    }

    updateCurrent(refs: StructureHierarchyRef[], action: 'add' | 'remove') {
        const hierarchy = this.current;
        const set = action === 'add'
            ? SetUtils.union(this.seletionSet, new Set(refs.map(r => r.cell.transform.ref)))
            : SetUtils.difference(this.seletionSet, new Set(refs.map(r => r.cell.transform.ref)));

        const trajectories = [];
        const models = [];
        const structures = [];

        for (const t of hierarchy.trajectories) {
            if (set.has(t.cell.transform.ref)) trajectories.push(t);
        }
        for (const m of hierarchy.models) {
            if (set.has(m.cell.transform.ref)) models.push(m);
        }
        for (const s of hierarchy.structures) {
            if (set.has(s.cell.transform.ref)) structures.push(s);
        }

        this._currentComponentGroups = void 0;
        this._currentSelectionSet = void 0;

        this.state.selection.trajectories = trajectories;
        this.state.selection.models = models;
        this.state.selection.structures = structures;

        this.behaviors.selection.next({ hierarchy, trajectories, models, structures });
    }

    remove(refs: (StructureHierarchyRef | string)[], canUndo?: boolean) {
        if (refs.length === 0) return;
        const deletes = this.plugin.state.data.build();
        for (const r of refs) deletes.delete(typeof r === 'string' ? r : r.cell.transform.ref);
        return deletes.commit({ canUndo: canUndo ? 'Remove' : false });
    }

    toggleVisibility(refs: ReadonlyArray<StructureHierarchyRef>, action?: 'show' | 'hide') {
        if (refs.length === 0) return;

        const isHidden = action !== void 0
            ? (action === 'show' ? false : true)
            : !refs[0].cell.state.isHidden;
        for (const c of refs) {
            setSubtreeVisibility(this.dataState, c.cell.transform.ref, isHidden);
        }
    }

    applyPreset<P = any, S = {}>(trajectories: ReadonlyArray<TrajectoryRef>, provider: TrajectoryHierarchyPresetProvider<P, S>, params?: P): Promise<any> {
        return this.plugin.dataTransaction(async () => {
            for (const t of trajectories) {
                if (t.models.length > 0) {
                    await this.clearTrajectory(t);
                }
                await this.plugin.builders.structure.hierarchy.applyPreset(t.cell, provider, params);
            }
        });
    }

    async updateStructure(s: StructureRef, params: any) {
        await this.plugin.dataTransaction(async () => {
            const root = StateTree.getDecoratorRoot(this.dataState.tree, s.cell.transform.ref);
            const children = this.dataState.tree.children.get(root).toArray();
            await this.remove(children, false);
            await this.plugin.state.updateTransform(this.plugin.state.data, s.cell.transform.ref, params, 'Structure Type');
            await this.plugin.builders.structure.representation.applyPreset(s.cell.transform.ref, 'auto');
        }, { canUndo: 'Structure Type' });
        PluginCommands.Camera.Reset(this.plugin);
    }

    private clearTrajectory(trajectory: TrajectoryRef) {
        const builder = this.dataState.build();
        for (const m of trajectory.models) {
            builder.delete(m.cell);
        }
        return builder.commit();
    }

    constructor(private plugin: PluginContext) {
        super();

        this.subscribe(plugin.state.data.events.changed, e => {
            if (e.inTransaction || plugin.behaviors.state.isAnimating.value) return;
            this.sync(true);
        });

        this.subscribe(plugin.behaviors.state.isAnimating, isAnimating => {
            if (!isAnimating && !plugin.behaviors.state.isUpdating.value) this.sync(true);
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

    export function getSelectedStructuresDescription(plugin: PluginContext) {
        const { structures } = plugin.managers.structure.hierarchy.selection;
        if (structures.length === 0) return '';

        if (structures.length === 1) {
            const s = structures[0];
            const data = s.cell.obj?.data;

            if (!data) return s.cell.obj?.label || 'Structure';

            const model = data.models[0] || data.representativeModel || data.masterModel;
            if (!model) return s.cell.obj?.label || 'Structure';

            const entryId = model.entryId;
            if (s.model?.trajectory?.models && s.model.trajectory.models.length === 1) return entryId;
            if (s.model) return `${s.model.cell.obj?.label} | ${entryId}`;
            return entryId;
        }

        const p = structures[0];
        const t = p?.model?.trajectory;
        let sameTraj = true;
        for (const s of structures) {
            if (s?.model?.trajectory !== t) {
                sameTraj = false;
                break;
            }
        }

        return sameTraj && t ? `${t.cell.obj?.label} | ${structures.length} structures` : `${structures.length} structures`;
    }
}