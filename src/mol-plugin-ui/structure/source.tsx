/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { HierarchyRef, ModelRef, TrajectoryRef } from '../../mol-plugin-state/manager/structure/hierarchy-state';
import { CollapsableControls, CollapsableState } from '../base';
import { ActionMenu } from '../controls/action-menu';
import { IconButton } from '../controls/common';
import { ParameterControls } from '../controls/parameters';
import { PluginCommands } from '../../mol-plugin/commands';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { memoize1 } from '../../mol-util/memoize';
import { StructureHierarchyManager } from '../../mol-plugin-state/manager/structure/hierarchy';

interface StructureSourceControlState extends CollapsableState {
    isBusy: boolean,
    show?: 'hierarchy' | 'actions'
}

export class StructureSourceControls extends CollapsableControls<{}, StructureSourceControlState> {
    protected defaultState(): StructureSourceControlState {
        return {
            header: 'Structure',
            isCollapsed: false,
            isBusy: false
        };
    }

    componentDidMount() {
        this.subscribe(this.plugin.managers.structure.hierarchy.behaviors.changed, () => this.forceUpdate());
        this.subscribe(this.plugin.behaviors.state.isBusy, v => {
            this.setState({ isBusy: v })
        });
    }

    private item = (ref: HierarchyRef) => {
        const selected = this.plugin.managers.structure.hierarchy.seletionSet;

        let label;
        switch (ref.kind) {
            case 'model': {
                const model = ref.cell.obj?.data;
                if (model?.trajectoryInfo.size! > 1) {
                    label = `${ref.cell.obj?.data.entryId} | Model ${model?.trajectoryInfo.index! + 1} of ${model?.trajectoryInfo.size}`;
                }
                label = `${ref.cell.obj?.data.entryId} | ${ref.cell.obj?.label}`; break;
            }
            case 'structure': {
                const model = ref.cell.obj?.data.models[0];
                if (model?.trajectoryInfo.size! > 1) {
                    label = `${ref.cell.obj?.data.models[0].entryId} | ${ref.cell.obj?.label} (Model ${model?.trajectoryInfo.index! + 1} of ${model?.trajectoryInfo.size})`; break;
                } else {
                    label = `${ref.cell.obj?.data.models[0].entryId} | ${ref.cell.obj?.label}`; break;
                }
            }
            default: label = ref.cell.obj?.label; break;
        }
        const item: ActionMenu.Item = { kind: 'item', label: label || ref.kind, selected: selected.has(ref.cell.transform.ref), value: [ref] };
        return item;
    }

    getTrajectoryItems = (t: TrajectoryRef): ActionMenu.Items => {
        if (t.models.length === 0) return this.item(t);
        return [ActionMenu.Header(t.cell.obj?.label!), ...t.models.map(this.getModelItems)];
    }

    private getModelItems = (m: ModelRef): ActionMenu.Items => {
        if (m.structures.length === 0) return this.item(m);
        if (m.structures.length === 1) {
            const selected = this.plugin.managers.structure.hierarchy.seletionSet;
            const ref = m.structures[0];
            return { label: `${m.cell.obj?.label} | ${ref.cell.obj?.label}`, selected: selected.has(ref.cell.transform.ref), value: [m, ref] } as ActionMenu.Item;
        }
        return [ActionMenu.Header(m.cell.obj?.label!), ...m.structures.map(this.item)];
    }

    get hierarchyItems() {
        const mng = this.plugin.managers.structure.hierarchy;
        const { current } = mng;
        const ret: ActionMenu.Items = [];

        if (current.trajectories.length > 1) {
            ret.push([
                ActionMenu.Header('Trajectories'),
                ...current.trajectories.map(this.item)
            ]);
        }

        if (current.models.length > 1) {
            ret.push([
                ActionMenu.Header('Models'),
                ...current.models.map(this.item)
            ])
        }

        if (current.trajectories.length === 1 && current.models.length === 1) {
            ret.push(...current.structures.map(this.item));
        } else if (current.structures.length > 0) {
            ret.push([
                ActionMenu.Header('Structures'),
                ...current.structures.map(this.item)
            ]);
        }

        return ret;
    }

    get isEmpty() {
        const { structures, models, trajectories } = this.plugin.managers.structure.hierarchy.current;
        return trajectories.length === 0 && models.length === 0 && structures.length === 0;
    }

    get label() {
        const { structures, models, trajectories } = this.plugin.managers.structure.hierarchy.state.selection;

        if (structures.length === 1) {
            const s = structures[0];
            if (s.model?.trajectory?.models && s.model.trajectory.models.length === 1) return s.cell.obj?.data.label;
            if (s.model) return `${s.model.cell.obj?.label} | ${s.cell.obj?.data.label}`;
            return s.cell.obj?.data.label;
        }

        if (structures.length > 1) {
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

        if (models.length > 0) {
            const t = models[0].trajectory;

            if (models.length === 1) {
                const model = models[0].cell.obj?.data;
                if (model?.trajectoryInfo.size! > 1) {
                    return `${t?.cell.obj?.label} | Model ${model?.trajectoryInfo.index! + 1} of ${model?.trajectoryInfo.size}`
                } else {
                    return `${t?.cell.obj?.label} | Model`
                }
            }

            let sameTraj = true;
            for (const m of models) {
                if (m.trajectory !== t) {
                    sameTraj = false;
                    break;
                }
            }

            return sameTraj ? `${t?.cell.obj?.label} | ${models.length} models` : `${models.length} models`;
        }

        if (trajectories.length > 0) {
            return trajectories.length === 1 ? `${trajectories[0].cell.obj?.label} trajectory` : `${trajectories.length} trajectories`;
        }

        if (trajectories.length === 0 && models.length === 0 && structures.length === 0) {
            return 'Nothing Loaded';
        }

        return 'Nothing Selected';
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

    actions = memoize1((sel: StructureHierarchyManager['selection']) => this._actions);

    get _actions() {
        const ret: ActionMenu.Items = [];

        const { selection } = this.plugin.managers.structure.hierarchy;
        if (selection.trajectories.some(t => t.cell.obj?.data.length! > 1)) {
            ret.push(ActionMenu.Item('Load all models', () => this.plugin.managers.structure.hierarchy.createModels(selection.trajectories, 'all')));
        }
        if (selection.trajectories.some(t => t.models.length > 1)) {
            ret.push(ActionMenu.Item('Load single model', () => this.plugin.managers.structure.hierarchy.createModels(selection.trajectories, 'single')));
        }
        
        // TODO: remove actions?
        return ret;
    }

    selectAction: ActionMenu.OnSelect = item => {
        if (!item) return;
        this.setState({ show: void 0 });
        (item?.value as any)();
    }

    updateStructureModel = async (params: any) => {
        const { selection } = this.plugin.managers.structure.hierarchy;
        const m = selection.structures[0].model!;
        this.plugin.state.updateTransform(this.plugin.state.data, m.cell.transform.ref, params, 'Model Index');
        // TODO: ?? PluginCommands.Camera.Reset(this.plugin);
    }

    get modelIndex() {
        const { selection } = this.plugin.managers.structure.hierarchy;
        if (selection.structures.length !== 1) return null;
        const m = selection.structures[0].model;
        if (!m || m.cell.transform.transformer !== StateTransforms.Model.ModelFromTrajectory) return null;
        if (m.cell.obj?.data.trajectoryInfo.size! <= 1) return null;

        const params = m.cell.params?.definition;
        if (!params) return null;

        return <ParameterControls params={params} values={m.cell.params?.values} onChangeObject={this.updateStructureModel} isDisabled={this.state.isBusy} />
    }

    updateStructure = async (params: any) => {
        const { selection } = this.plugin.managers.structure.hierarchy;
        const s = selection.structures[0];
        await this.plugin.state.updateTransform(this.plugin.state.data, s.cell.transform.ref, params, 'Structure Type');
        PluginCommands.Camera.Reset(this.plugin);
    }

    get structureType() {
        const { selection } = this.plugin.managers.structure.hierarchy;
        if (selection.structures.length !== 1) return null;
        const s = selection.structures[0];
        const params = s.cell.params?.definition;
        if (!params) return null;

        return <ParameterControls params={params} values={s.cell.params?.values} onChangeObject={this.updateStructure} isDisabled={this.state.isBusy} />
    }

    renderControls() {
        const disabled = this.state.isBusy || this.isEmpty;
        const actions = this.actions(this.plugin.managers.structure.hierarchy.selection);
        return <>
            <div className='msp-btn-row-group' style={{ marginTop: '1px' }}>
                <button className='msp-btn msp-form-control msp-flex-item' onClick={this.toggleHierarchy} style={{ overflow: 'hidden', textOverflow: 'ellipsis' }} disabled={disabled}>
                    {this.label}
                </button>
                {actions.length > 0 && <IconButton customClass='msp-form-control' style={{ flex: '0 0 32px' }} onClick={this.toggleActions} icon='dot-3' title='Actions' toggleState={this.state.show === 'actions'} disabled={disabled} />}
            </div>
            {this.state.show === 'hierarchy' && <ActionMenu items={this.hierarchyItems} onSelect={this.selectHierarchy} multiselect />}
            {this.state.show === 'actions' && <ActionMenu items={actions} onSelect={this.selectAction} />}
            {this.modelIndex}
            {this.structureType}
        </>;
    }
}