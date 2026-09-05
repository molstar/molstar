/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { setSubtreeVisibility } from '../../../mol-plugin/behavior/static/state';
import { PluginContext } from '../../../mol-plugin/context';
import { PluginComponent } from '../../component';
import { buildParticleHierarchy, ParticleHierarchy, ParticleHierarchyRef, ParticleListRef } from './hierarchy-state';
import { StateTransforms } from '../../transforms';

export class ParticleHierarchyManager extends PluginComponent {
    private state = {
        syncedTree: this.dataState.tree,
        notified: false,

        hierarchy: ParticleHierarchy(),
        selection: void 0 as ParticleListRef | undefined
    };

    readonly behaviors = {
        selection: this.ev.behavior({
            hierarchy: this.current,
            list: this.selection
        })
    };

    private get dataState() {
        return this.plugin.state.data;
    }

    get current() {
        this.sync(false);
        return this.state.hierarchy;
    }

    get selection() {
        this.sync(false);
        return this.state.selection;
    }

    private sync(notify: boolean) {
        if (!notify && this.dataState.inUpdate) return;

        if (this.state.syncedTree === this.dataState.tree) {
            if (notify && !this.state.notified) {
                this.state.notified = true;
                this.behaviors.selection.next({ hierarchy: this.state.hierarchy, list: this.state.selection });
            }

            return;
        }

        this.state.syncedTree = this.dataState.tree;

        const update = buildParticleHierarchy(this.plugin.state.data, this.current);
        if (!update.changed) {
            return;
        }

        const { hierarchy } = update;

        this.state.hierarchy = hierarchy;
        if (!this.state.selection) {
            this.state.selection = hierarchy.lists[0];
        } else {
            this.state.selection = hierarchy.refs.has(this.state.selection.cell.transform.ref) ? hierarchy.refs.get(this.state.selection.cell.transform.ref) as ParticleListRef : hierarchy.lists[0];
        }

        if (notify) {
            this.state.notified = true;
            this.behaviors.selection.next({ hierarchy, list: this.state.selection });
        } else {
            this.state.notified = false;
        }
    }

    setCurrent(list?: ParticleListRef) {
        this.state.selection = list || this.state.hierarchy.lists[0];
        this.behaviors.selection.next({ hierarchy: this.state.hierarchy, list: list || this.state.hierarchy.lists[0] });
    }

    remove(refs: (ParticleHierarchyRef | string)[], canUndo?: boolean) {
        if (refs.length === 0) return;
        const deletes = this.plugin.state.data.build();
        for (const r of refs) deletes.delete(typeof r === 'string' ? r : r.cell.transform.ref);
        return deletes.commit({ canUndo: canUndo ? 'Remove' : false });
    }

    toggleVisibility(refs: ReadonlyArray<ParticleHierarchyRef>, action?: 'show' | 'hide') {
        if (refs.length === 0) return;

        const isHidden = action !== void 0
            ? (action === 'show' ? false : true)
            : !refs[0].cell.state.isHidden;
        for (const c of refs) {
            setSubtreeVisibility(this.dataState, c.cell.transform.ref, isHidden);
        }
    }

    addRepresentation(ref: ParticleListRef, type: string) {
        const update = this.dataState.build()
            .to(ref.cell)
            .apply(StateTransforms.Particles.ParticlesRepresentation3D, { type: { name: type, params: {} } });

        return update.commit({ canUndo: 'Add Representation' });
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

export namespace ParticleHierarchyManager {
    export function getRepresentationTypes(plugin: PluginContext, pivot: ParticleListRef | undefined) {
        return pivot?.cell.obj?.data
            ? plugin.representation.particles.registry.getApplicableTypes(pivot.cell.obj?.data!)
            : plugin.representation.particles.registry.types;
    }
}
