/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { setSubtreeVisibility } from '../../../mol-plugin/behavior/static/state';
import { PluginContext } from '../../../mol-plugin/context';
import { PluginComponent } from '../../component';
import { buildVolumeHierarchy, VolumeHierarchy, VolumeHierarchyRef, VolumeRef } from './hierarchy-state';
import { createVolumeRepresentationParams } from '../../helpers/volume-representation-params';
import { StateTransforms } from '../../transforms';

export class VolumeHierarchyManager extends PluginComponent {
    private state = {
        syncedTree: this.dataState.tree,
        notified: false,

        hierarchy: VolumeHierarchy(),
        selection: void 0 as VolumeRef | undefined
    }

    readonly behaviors = {
        selection: this.ev.behavior({
            hierarchy: this.current,
            volume: this.selection
        })
    }

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
                this.behaviors.selection.next({ hierarchy: this.state.hierarchy, volume: this.state.selection });
            }

            return;
        }

        this.state.syncedTree = this.dataState.tree;

        const update = buildVolumeHierarchy(this.plugin.state.data, this.current);
        if (!update.changed) {
            return;
        }

        const { hierarchy } = update;

        this.state.hierarchy = hierarchy;
        if (!this.state.selection) {
            this.state.selection = hierarchy.volumes[0];
        } else {
            this.state.selection = hierarchy.refs.has(this.state.selection.cell.transform.ref) ? hierarchy.refs.get(this.state.selection.cell.transform.ref) as VolumeRef : hierarchy.volumes[0];
        }

        if (notify) {
            this.state.notified = true;
            this.behaviors.selection.next({ hierarchy, volume: this.state.selection });
        } else {
            this.state.notified = false;
        }
    }

    setCurrent(volume?: VolumeRef) {
        this.state.selection = volume || this.state.hierarchy.volumes[0];
        this.behaviors.selection.next({ hierarchy: this.state.hierarchy, volume: volume || this.state.hierarchy.volumes[0] });
    }

    // TODO: have common util
    remove(refs: (VolumeHierarchyRef | string)[], canUndo?: boolean) {
        if (refs.length === 0) return;
        const deletes = this.plugin.state.data.build();
        for (const r of refs) deletes.delete(typeof r === 'string' ? r : r.cell.transform.ref);
        return deletes.commit({ canUndo: canUndo ? 'Remove' : false });
    }

    // TODO: have common util
    toggleVisibility(refs: ReadonlyArray<VolumeHierarchyRef>, action?: 'show' | 'hide') {
        if (refs.length === 0) return;

        const isHidden = action !== void 0
            ? (action === 'show' ? false : true)
            : !refs[0].cell.state.isHidden;
        for (const c of refs) {
            setSubtreeVisibility(this.dataState, c.cell.transform.ref, isHidden);
        }
    }

    addRepresentation(ref: VolumeRef, type: string) {
        const update = this.dataState.build()
            .to(ref.cell)
            .apply(StateTransforms.Representation.VolumeRepresentation3D, createVolumeRepresentationParams(this.plugin, ref.cell.obj?.data, {
                type: type as any,
            }));

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

export namespace VolumeHierarchyManager {
    export function getRepresentationTypes(plugin: PluginContext, pivot: VolumeRef | undefined) {
        return pivot?.cell.obj?.data
            ? plugin.representation.volume.registry.getApplicableTypes(pivot.cell.obj?.data!)
            : plugin.representation.volume.registry.types;
    }
}