/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '../../mol-plugin/context';
import { StructureHierarchy, buildStructureHierarchy, ModelRef, StructureComponentRef } from './structure/hierarchy';
import { RxEventHelper } from '../../mol-util/rx-event-helper';

export class StructureHierarchyManager {
    private ev = RxEventHelper.create();

    readonly behaviors = {
        hierarchy: this.ev.behavior(StructureHierarchy()),
        currentModels: this.ev.behavior([] as ReadonlyArray<ModelRef>)
    }

    private checkCurrent() {
        const hierarchy = this.behaviors.hierarchy.value;
        const current = this.behaviors.currentModels.value;
        if (current.length === 0) {
            const models = hierarchy.trajectories[0]?.models;
            if (models) {
                this.behaviors.currentModels.next(models);
            }
            return;
        }

        const newCurrent: ModelRef[] = [];
        for (const c of current) {
            const ref = hierarchy.refs.get(c.cell.transform.ref) as ModelRef;
            if (!ref) continue;
            newCurrent.push(ref);
        }

        if (newCurrent.length === 0 && hierarchy.trajectories[0]?.models) {
            this.behaviors.currentModels.next(hierarchy.trajectories[0]?.models);
        }

        this.behaviors.currentModels.next(newCurrent);
    }

    private sync() {
        const update = buildStructureHierarchy(this.plugin.state.dataState, this.behaviors.hierarchy.value);
        if (update.added.length === 0 && update.updated.length === 0 && update.removed.length === 0) {
            return;
        }

        this.behaviors.hierarchy.next(update.hierarchy)
        this.checkCurrent();
    }

    constructor(private plugin: PluginContext) {
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
    function componentKey(c: StructureComponentRef) {
        if (!c.cell.transform.tags) return;
        return [...c.cell.transform.tags].sort().join();
    }

    export function getCommonComponentPivots(models: ReadonlyArray<ModelRef>) {
        if (!models[0]?.structures?.length) return [];
        if (models[0]?.structures?.length === 1) return models[0]?.structures[0]?.components || [];

        const pivots = new Map<string, StructureComponentRef>();

        for (const c of models[0]?.structures[0]?.components) {
            const key = componentKey(c);
            if (!key) continue;
            pivots.set(key, c);
        }

        for (const m of models) {
            for (const s of m.structures) {
                for (const c of s.components) {
                    const key = componentKey(c);
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