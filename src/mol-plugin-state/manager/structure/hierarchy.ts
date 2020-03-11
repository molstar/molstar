/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '../../../mol-plugin/context';
import { StructureHierarchy, buildStructureHierarchy, ModelRef, StructureComponentRef } from './hierarchy-state';
import { PluginComponent } from '../../component';

interface StructureHierarchyManagerState {
    hierarchy: StructureHierarchy,
    currentModels: ReadonlyArray<ModelRef>,
}

export class StructureHierarchyManager extends PluginComponent<StructureHierarchyManagerState> {
    readonly behaviors = {
        hierarchy: this.ev.behavior(this.state.hierarchy),
        currentModels: this.ev.behavior(this.state.currentModels)
    }

    private syncCurrent(hierarchy: StructureHierarchy) {
        const current = this.behaviors.currentModels.value;
        if (current.length === 0) {
            const models = hierarchy.trajectories[0]?.models;
            if (models) {
                return models;
            }
            return [];
        }

        const newCurrent: ModelRef[] = [];
        for (const c of current) {
            const ref = hierarchy.refs.get(c.cell.transform.ref) as ModelRef;
            if (!ref) continue;
            newCurrent.push(ref);
        }

        if (newCurrent.length === 0 && hierarchy.trajectories[0]?.models) {
            return hierarchy.trajectories[0]?.models;
        }

        return newCurrent;
    }

    private sync() {
        const update = buildStructureHierarchy(this.plugin.state.dataState, this.behaviors.hierarchy.value);
        if (update.added.length === 0 && update.updated.length === 0 && update.removed.length === 0) {
            return;
        }

        const currentModels = this.syncCurrent(update.hierarchy);
        this.updateState({ hierarchy: update.hierarchy, currentModels });

        this.behaviors.hierarchy.next(this.state.hierarchy);
        this.behaviors.currentModels.next(this.state.currentModels);
    }

    constructor(private plugin: PluginContext) {
        super({
            hierarchy: StructureHierarchy(),
            currentModels: []
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
    export function getCommonComponentPivots(models: ReadonlyArray<ModelRef>) {
        if (!models[0]?.structures?.length) return [];
        if (models[0]?.structures?.length === 1) return models[0]?.structures[0]?.components || [];

        const pivots = new Map<string, StructureComponentRef>();

        for (const c of models[0]?.structures[0]?.components) {
            const key = c.key;
            if (!key) continue;
            pivots.set(key, c);
        }

        for (const m of models) {
            for (const s of m.structures) {
                for (const c of s.components) {
                    const key = c.key;
                    if (!key) continue;
                    if (!pivots.has(key)) pivots.delete(key);
                }
            }
        }

        const ret: StructureComponentRef[] = [];
        pivots.forEach(function (this: StructureComponentRef[], p) { this.push(p) }, ret);
        return ret;
    }
}