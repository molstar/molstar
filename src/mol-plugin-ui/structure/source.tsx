/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { StructureHierarchyRef, ModelRef, TrajectoryRef } from '../../mol-plugin-state/manager/structure/hierarchy-state';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { CollapsableControls, CollapsableState } from '../base';
import { ActionMenu } from '../controls/action-menu';
import { Button, IconButton, ExpandGroup } from '../controls/common';
import { ParameterControls } from '../controls/parameters';
import { StructureFocusControls } from './focus';
import { UpdateTransformControl } from '../state/update-transform';
import { StructureSelectionStatsControls } from './selection';
import { StateSelection } from '../../mol-state';
import { MoleculeSvg, BookmarksOutlinedSvg } from '../controls/icons';
import { Model } from '../../mol-model/structure';

interface StructureSourceControlState extends CollapsableState {
    isBusy: boolean,
    show?: 'hierarchy' | 'presets'
}

export class StructureSourceControls extends CollapsableControls<{}, StructureSourceControlState> {
    protected defaultState(): StructureSourceControlState {
        return {
            header: 'Structure',
            isCollapsed: false,
            isBusy: false,
            brand: { accent: 'purple', svg: MoleculeSvg }
        };
    }

    componentDidMount() {
        this.subscribe(this.plugin.managers.structure.hierarchy.behaviors.selection, () => this.forceUpdate());
        this.subscribe(this.plugin.behaviors.state.isBusy, v => {
            this.setState({ isBusy: v });
        });
    }

    private item = (ref: StructureHierarchyRef) => {
        const selected = this.plugin.managers.structure.hierarchy.seletionSet;

        let label;
        switch (ref.kind) {
            case 'model': {
                const model = ref.cell.obj?.data;
                if (model && Model.TrajectoryInfo.get(model).size > 1) {
                    label = `${ref.cell.obj?.data.entryId} | Model ${Model.TrajectoryInfo.get(model).index + 1} of ${Model.TrajectoryInfo.get(model).size}`;
                }
                label = `${ref.cell.obj?.data.entryId} | ${ref.cell.obj?.label}`; break;
            }
            case 'structure': {
                const model = ref.cell.obj?.data.models[0];
                if (model && Model.TrajectoryInfo.get(model).size! > 1) {
                    label = `${model.entryId} | ${ref.cell.obj?.label} (Model ${Model.TrajectoryInfo.get(model).index + 1} of ${Model.TrajectoryInfo.get(model).size})`; break;
                } else if (model) {
                    label = `${model.entryId} | ${ref.cell.obj?.label}`; break;
                } else {
                    label = `${ref.cell.obj?.label}`; break;
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

        if (current.models.length > 1 || current.trajectories.length > 1) {
            ret.push([
                ActionMenu.Header('Models'),
                ...current.models.map(this.item)
            ]);
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
        const { structures, models, trajectories } = this.plugin.managers.structure.hierarchy.selection;

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
                if (model && Model.TrajectoryInfo.get(model).size > 1) {
                    return `${t?.cell.obj?.label} | Model ${Model.TrajectoryInfo.get(model).index + 1} of ${Model.TrajectoryInfo.get(model).size}`;
                } else {
                    return `${t?.cell.obj?.label} | Model`;
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
        if (!items || items.length === 0) return;

        const refs: StructureHierarchyRef[] = [];
        for (const i of items) {
            for (const r of (i.value as StructureHierarchyRef[])) refs.push(r);
        }

        this.plugin.managers.structure.hierarchy.updateCurrent(refs, items[0].selected ? 'remove' : 'add');
    }

    toggleHierarchy = () => this.setState({ show: this.state.show !== 'hierarchy' ? 'hierarchy' : void 0 });
    togglePreset = () => this.setState({ show: this.state.show !== 'presets' ? 'presets' : void 0 });

    get presetActions() {
        const actions: ActionMenu.Item[] = [];
        const { trajectories } = this.plugin.managers.structure.hierarchy.selection;
        if (trajectories.length !== 1) return actions;

        const providers = this.plugin.builders.structure.hierarchy.getPresets(trajectories[0].cell.obj);
        for (const p of providers) {
            actions.push(ActionMenu.Item(p.display.name, p, { description: p.display.description }));
        }
        return actions;
    }

    applyPreset: ActionMenu.OnSelect = item => {
        this.setState({ show: void 0 });

        if (!item) return;
        const mng = this.plugin.managers.structure;

        const { trajectories } = mng.hierarchy.selection;
        mng.hierarchy.applyPreset(trajectories, item.value as any);
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
        if (!m.cell.obj || Model.TrajectoryInfo.get(m.cell.obj.data).size <= 1) return null;

        const params = m.cell.params?.definition;
        if (!params) return null;

        return <ParameterControls params={params} values={m.cell.params?.values} onChangeValues={this.updateStructureModel} isDisabled={this.state.isBusy} />;
    }

    updateStructure = (params: any) => {
        const { selection } = this.plugin.managers.structure.hierarchy;
        const s = selection.structures[0];
        return this.plugin.managers.structure.hierarchy.updateStructure(s, params);
    }

    get structureType() {
        const { selection } = this.plugin.managers.structure.hierarchy;
        if (selection.structures.length !== 1) return null;
        const s = selection.structures[0];
        const params = s.cell.params?.definition;
        if (!params || !s.cell.parent) return null;

        return <UpdateTransformControl state={s.cell.parent} transform={s.cell.transform} customHeader='none' customUpdate={this.updateStructure} noMargin autoHideApply />;
    }

    get transform() {
        const { selection } = this.plugin.managers.structure.hierarchy;
        if (selection.structures.length !== 1) return null;
        const pivot = selection.structures[0];
        if (!pivot.cell.parent) return null;
        const t = StateSelection.tryFindDecorator(this.plugin.state.data, pivot.cell.transform.ref, StateTransforms.Model.TransformStructureConformation);
        if (!t) return;

        return <ExpandGroup header={`Conformation Transform`}>
            <UpdateTransformControl state={t.parent!} transform={t.transform} customHeader='none' noMargin autoHideApply />
        </ExpandGroup>;
    }

    renderControls() {
        const disabled = this.state.isBusy || this.isEmpty;
        const presets = this.presetActions;
        const label = this.label;
        return <>
            <div className='msp-flex-row' style={{ marginTop: '1px' }}>
                <Button noOverflow flex onClick={this.toggleHierarchy} disabled={disabled} title={label}>{label}</Button>
                {presets.length > 0 && <IconButton svg={BookmarksOutlinedSvg} className='msp-form-control' flex='40px' onClick={this.togglePreset} title='Apply a structure presets to the current hierarchy.' toggleState={this.state.show === 'presets'} disabled={disabled} />}
            </div>
            {this.state.show === 'hierarchy' && <ActionMenu items={this.hierarchyItems} onSelect={this.selectHierarchy} multiselect />}
            {this.state.show === 'presets' && <ActionMenu items={presets} onSelect={this.applyPreset} />}
            {this.modelIndex}
            {this.structureType}
            {this.transform}

            <div style={{ marginTop: '6px' }}>
                <StructureFocusControls />
                <StructureSelectionStatsControls hideOnEmpty />
            </div>
        </>;
    }
}