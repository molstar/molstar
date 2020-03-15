/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { HierarchyRef, ModelRef, TrajectoryRef } from '../../mol-plugin-state/manager/structure/hierarchy-state';
import { CollapsableControls, CollapsableState } from '../base';
import { ActionMenu } from '../controls/action-menu';
import { Icon } from '../controls/icons';

interface StructureSourceControlState extends CollapsableState {
    isBusy: boolean,
    show?: 'hierarchy' | 'actions'
}

export class StructureSourceControls extends CollapsableControls<{}, StructureSourceControlState> {
    protected defaultState(): StructureSourceControlState {
        return { 
            header: 'Source',
            isCollapsed: false,
            isBusy: false
        };
    }

    componentDidMount() {
        this.subscribe(this.plugin.managers.structure.hierarchy.behaviors.current, () => this.forceUpdate());
        this.subscribe(this.plugin.behaviors.state.isBusy, v => {
            this.setState({ isBusy: v })
        });
    }

    private item = (ref: HierarchyRef) => {
        const selected = this.plugin.managers.structure.hierarchy.currentSeletionSet;
        return { label: ref.cell.obj?.label, selected: selected.has(ref.cell.transform.ref), value: [ref] } as ActionMenu.Item;
    }

    getTrajectoryItems = (t: TrajectoryRef): ActionMenu.Items => {
        if (t.models.length === 0) return this.item(t);
        return [ActionMenu.Header(t.cell.obj?.label!), ...t.models.map(this.getModelItems)];
    }

    private getModelItems = (m: ModelRef): ActionMenu.Items => {
        if (m.structures.length === 0) return this.item(m);
        if (m.structures.length === 1) {
            const selected = this.plugin.managers.structure.hierarchy.currentSeletionSet;
            const ref = m.structures[0];
            return { label: `${m.cell.obj?.label} | ${ref.cell.obj?.label}`, selected: selected.has(ref.cell.transform.ref), value: [m, ref] } as ActionMenu.Item;
        }
        return [ActionMenu.Header(m.cell.obj?.label!), ...m.structures.map(this.item)];
    }

    get hierarchyItems() {
        return this.plugin.managers.structure.hierarchy.state.current.trajectories.map(this.getTrajectoryItems);
    }

    get label() {
        const { structures, models, trajectories } = this.plugin.managers.structure.hierarchy.state.current;

        // TODO: better labels
        
        if (structures.length === 1) {
            const s = structures[0];
            if (s.model?.trajectory?.models && s.model.trajectory.models.length === 1) return s.cell.obj?.data.label;
            if (s.model) return `${s.model.cell.obj?.label} | ${s.cell.obj?.data.label}`;
            return s.cell.obj?.data.label;
        }

        if (structures.length > 1) {
            return `${structures.length} structures`;
        }

        if (models.length > 0) {
            return `${models.length} model(s)`;
        }

        if (trajectories.length > 0) {
            return `${trajectories.length} trajectory(s)`;
        }

        return 'No structure loaded';
    }

    selectHierarchy: ActionMenu.OnSelectMany = (items) => {
        if (!items || items.length === 0) return 0;

        const refs: HierarchyRef[] = [];
        for (const i of items) {
            for (const r of (i.value as HierarchyRef[])) refs.push(r);
        }

        this.plugin.managers.structure.hierarchy.updateCurrent(refs, items[0].selected ? 'remove' : 'add')
    }

    toggleHierarchy = () => this.setState({ show: this.state.show !== 'hierarchy' ? 'hierarchy' : void 0 });
    toggleActions = () => this.setState({ show: this.state.show !== 'actions' ? 'actions' : void 0 });

    get actions() {
        const ret: ActionMenu.Items = [
            ActionMenu.Item('Show all models', () => this.plugin.managers.structure.hierarchy.createAllModels(this.plugin.managers.structure.hierarchy.state.current.trajectories[0]))
        ];
        return ret;
    }

    selectAction: ActionMenu.OnSelect = item => {
        if (!item) return;
        this.setState({ show: void 0 });
        (item?.value as any)();
    }

    renderControls() {
        return <>
            <div className='msp-flex-row' style={{ marginTop: '1px' }}>
                <button className='msp-btn msp-form-control msp-flex-item' onClick={this.toggleHierarchy} style={{ overflow: 'hidden', textOverflow: 'ellipsis' }}>
                    {this.label}
                </button>
                <button className='msp-btn msp-form-control msp-flex-item' onClick={this.toggleActions} style={{ flex: '0 0 40px', }}>
                    <Icon name='dot-3' />
                </button>
            </div>
            {this.state.show === 'hierarchy' && <ActionMenu items={this.hierarchyItems} onSelect={this.selectHierarchy} multiselect />}
            {this.state.show === 'actions' && <ActionMenu items={this.actions} onSelect={this.selectAction} />}
        </>;
    }
}