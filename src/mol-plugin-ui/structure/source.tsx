/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { HierarchyRef, ModelRef, TrajectoryRef } from '../../mol-plugin-state/manager/structure/hierarchy-state';
import { CollapsableControls, CollapsableState } from '../base';
import { ActionMenu } from '../controls/action-menu';

interface StructureSourceControlState extends CollapsableState {
    isBusy: boolean
}

export class StructureSourceControls extends CollapsableControls<{}, StructureSourceControlState> {
    protected defaultState(): StructureSourceControlState {
        return { header: 'Source', isCollapsed: false, isBusy: false };
    }

    componentDidMount() {
        this.subscribe(this.plugin.managers.structure.hierarchy.behaviors.current, () => this.forceUpdate());
        this.subscribe(this.plugin.behaviors.state.isBusy, v => {
            this.setState({ isBusy: v })
        });
    }

    private item = (ref: HierarchyRef) => {
        const selected = this.plugin.managers.structure.hierarchy.currentSeletionSet;
        return { label: ref.cell.obj?.label, selected: selected.has(ref.cell.transform.ref), value: ref } as ActionMenu.Item;
    }

    getTrajectoryItems = (t: TrajectoryRef): ActionMenu.Items => {
        if (t.models.length === 0) return this.item(t);
        // if (t.models.length === 1) return this.getModelItems(t.models[0]);
        return [t.cell.obj?.label!, ...t.models.map(this.getModelItems)];
    }

    private getModelItems = (m: ModelRef): ActionMenu.Items => {
        if (m.structures.length === 0) return this.item(m);
        // if (m.structures.length === 1) return this.item(m.structures[0]);
        return [m.cell.obj?.label!, ...m.structures.map(this.item)];
    }

    get hierarchyItems() {
        return this.plugin.managers.structure.hierarchy.state.current.trajectories.map(this.getTrajectoryItems);
    }

    onSelect: ActionMenu.OnSelectMany = (items) => {
        if (!items || items.length === 0) return 0;
        this.plugin.managers.structure.hierarchy.updateCurrent(items.map(i => i.value as HierarchyRef), items[0].selected ? 'remove' : 'add')
    }

    renderControls() {
        return <>
            <ActionMenu items={this.hierarchyItems} onSelect={this.onSelect} multiselect />
            <button onClick={() => this.plugin.managers.structure.hierarchy.createAllModels(this.plugin.managers.structure.hierarchy.state.current.trajectories[0])}>All Models</button>
        </>;
    }
}