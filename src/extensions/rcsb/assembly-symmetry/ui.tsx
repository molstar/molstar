/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { CollapsableState, CollapsableControls } from '../../../mol-plugin-ui/base';
import { ApplyActionControl } from '../../../mol-plugin-ui/state/apply-action';
import { InitAssemblySymmetry3D, AssemblySymmetry3D, AssemblySymmetryPreset, tryCreateAssemblySymmetry } from './behavior';
import { AssemblySymmetryProvider,  AssemblySymmetryProps, AssemblySymmetryDataProvider, AssemblySymmetry } from './prop';
import { ParameterControls } from '../../../mol-plugin-ui/controls/parameters';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { StructureHierarchyManager } from '../../../mol-plugin-state/manager/structure/hierarchy';
import { StateAction, StateSelection } from '../../../mol-state';
import { PluginStateObject } from '../../../mol-plugin-state/objects';
import { PluginContext } from '../../../mol-plugin/context';
import { Task } from '../../../mol-task';
import { ExtensionSvg, CheckSvg } from '../../../mol-plugin-ui/controls/icons';

interface AssemblySymmetryControlState extends CollapsableState {
    isBusy: boolean
}

export class AssemblySymmetryControls extends CollapsableControls<{}, AssemblySymmetryControlState> {
    protected defaultState(): AssemblySymmetryControlState {
        return {
            header: 'Assembly Symmetry',
            isCollapsed: false,
            isBusy: false,
            isHidden: true,
            brand: { accent: 'cyan', svg: ExtensionSvg }
        };
    }

    componentDidMount() {
        this.subscribe(this.plugin.managers.structure.hierarchy.behaviors.selection, () => {
            this.setState({
                isHidden: !this.canEnable(),
                description: StructureHierarchyManager.getSelectedStructuresDescription(this.plugin)
            });
        });
        this.subscribe(this.plugin.state.events.cell.stateUpdated, e => {
            if (e.cell.transform.transformer === AssemblySymmetry3D) this.forceUpdate();
        });
        this.subscribe(this.plugin.behaviors.state.isBusy, v => this.setState({ isBusy: v }));
    }

    get pivot() {
        return this.plugin.managers.structure.hierarchy.selection.structures[0];
    }

    canEnable() {
        const { selection } = this.plugin.managers.structure.hierarchy;
        if (selection.structures.length !== 1) return false;
        const pivot = this.pivot.cell;
        if (!pivot.obj) return false;
        return !!InitAssemblySymmetry3D.definition.isApplicable?.(pivot.obj, pivot.transform, this.plugin);
    }

    renderEnable() {
        const pivot = this.pivot;
        if (!pivot.cell.parent) return null;
        return <ApplyActionControl state={pivot.cell.parent} action={EnableAssemblySymmetry3D} initiallyCollapsed={true} nodeRef={pivot.cell.transform.ref} simpleApply={{ header: 'Enable', icon: CheckSvg }} />;
    }

    renderNoSymmetries() {
        return <div className='msp-row-text'>
            <div>No Symmetries for Assembly</div>
        </div>;
    }

    get params() {
        const structure = this.pivot.cell.obj?.data;
        const params = PD.clone(structure ? AssemblySymmetryProvider.getParams(structure) : AssemblySymmetryProvider.defaultParams);
        params.serverUrl.isHidden = true;
        params.symmetryIndex.options = [[-1, 'Off'], ...params.symmetryIndex.options];
        return params;
    }

    get values() {
        const structure = this.pivot.cell.obj?.data;
        if (structure) {
            return AssemblySymmetryProvider.props(structure);
        } else {
            return { ...PD.getDefaultValues(AssemblySymmetryProvider.defaultParams), symmetryIndex: -1 };
        }
    }

    async updateAssemblySymmetry(values: AssemblySymmetryProps) {
        const s = this.pivot;
        const currValues = AssemblySymmetryProvider.props(s.cell.obj!.data);
        if (PD.areEqual(AssemblySymmetryProvider.defaultParams, currValues, values)) return;

        if (s.properties) {
            const b = this.plugin.state.data.build();
            b.to(s.properties.cell).update(old => {
                old.properties[AssemblySymmetryProvider.descriptor.name] = values;
            });
            await b.commit();
        } else {
            const pd = this.plugin.customStructureProperties.getParams(s.cell.obj?.data);
            const params = PD.getDefaultValues(pd);
            params.properties[AssemblySymmetryProvider.descriptor.name] = values;
            await this.plugin.builders.structure.insertStructureProperties(s.cell, params);
        }

        for (const components of this.plugin.managers.structure.hierarchy.currentComponentGroups) {
            if (values.symmetryIndex === -1) {
                const name = components[0]?.representations[0]?.cell.transform.params?.colorTheme.name;
                if (name === AssemblySymmetry.Tag.Cluster) {
                    await this.plugin.managers.structure.component.updateRepresentationsTheme(components, { color: 'default' });
                }
            } else {
                tryCreateAssemblySymmetry(this.plugin, s.cell);
                await this.plugin.managers.structure.component.updateRepresentationsTheme(components, { color: AssemblySymmetry.Tag.Cluster as any });
            }
        }
    }

    paramsOnChange = (options: AssemblySymmetryProps) => {
        this.updateAssemblySymmetry(options);
    }

    get hasAssemblySymmetry3D() {
        return !this.pivot.cell.parent || !!StateSelection.findTagInSubtree(this.pivot.cell.parent.tree, this.pivot.cell.transform.ref, AssemblySymmetry.Tag.Representation);
    }

    get enable() {
        return !this.hasAssemblySymmetry3D && this.values.symmetryIndex !== -1;
    }

    get noSymmetries() {
        const structure = this.pivot.cell.obj?.data;
        const data = structure && AssemblySymmetryDataProvider.get(structure).value;
        return data && data.filter(sym => sym.symbol !== 'C1').length === 0;
    }

    renderParams() {
        return <>
            <ParameterControls params={this.params} values={this.values} onChangeValues={this.paramsOnChange} />
        </>;
    }

    renderControls() {
        if (!this.pivot) return null;
        if (this.noSymmetries) return this.renderNoSymmetries();
        if (this.enable) return this.renderEnable();
        return this.renderParams();
    }
}

const EnableAssemblySymmetry3D = StateAction.build({
    from: PluginStateObject.Molecule.Structure,
})(({ a, ref, state }, plugin: PluginContext) => Task.create('Enable Assembly Symmetry', async ctx => {
    await AssemblySymmetryPreset.apply(ref, Object.create(null), plugin);
}));