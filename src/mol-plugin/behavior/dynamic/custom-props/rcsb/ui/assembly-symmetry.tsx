/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { CollapsableState, CollapsableControls } from '../../../../../../mol-plugin-ui/base';
import { ApplyActionControl } from '../../../../../../mol-plugin-ui/state/apply-action';
import { InitAssemblySymmetry3D, AssemblySymmetry3D } from '../assembly-symmetry';
import { AssemblySymmetryProvider,  AssemblySymmetryProps, AssemblySymmetryDataProvider } from '../../../../../../mol-model-props/rcsb/assembly-symmetry';
import { ParameterControls } from '../../../../../../mol-plugin-ui/controls/parameters';
import { ParamDefinition as PD } from '../../../../../../mol-util/param-definition';
import { StructureHierarchyManager } from '../../../../../../mol-plugin-state/manager/structure/hierarchy';

interface AssemblySymmetryControlState extends CollapsableState {
    isBusy: boolean
    options: AssemblySymmetryProps
}

export class AssemblySymmetryControls extends CollapsableControls<{}, AssemblySymmetryControlState> {
    protected defaultState(): AssemblySymmetryControlState {
        return {
            header: 'Assembly Symmetry',
            isCollapsed: false,
            isBusy: false,
            isHidden: true,
            options: PD.getDefaultValues(AssemblySymmetryProvider.defaultParams)
        };
    }

    componentDidMount() {
        this.subscribe(this.plugin.managers.structure.hierarchy.behaviors.selection, () => this.setState({ isHidden: !this.canEnable() }));
        this.subscribe(this.plugin.behaviors.state.isBusy, v => this.setState({ isBusy: v }));

        this.subscribe(this.plugin.managers.structure.hierarchy.behaviors.selection, c => this.setState({
            description: StructureHierarchyManager.getSelectedStructuresDescription(this.plugin)
        }));
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
        return <ApplyActionControl state={pivot.cell.parent} action={InitAssemblySymmetry3D} initiallyCollapsed={true} nodeRef={pivot.cell.transform.ref} simpleApply={{ header: 'Enable', icon: 'check' }} />;
    }

    renderNoSymmetries() {
        return <div className='msp-control-row msp-row-text'>
            <div>No Symmetries for Assembly</div>
        </div>;
    }

    get params() {
        const structure = this.pivot.cell.obj?.data
        const params = PD.clone(structure ? AssemblySymmetryProvider.getParams(structure) : AssemblySymmetryProvider.defaultParams)
        params.serverUrl.isHidden = true
        params.symmetryIndex.options = [[-1, 'Off'], ...params.symmetryIndex.options]
        return params
    }

    async updateAssemblySymmetry(values: AssemblySymmetryProps) {
        const s = this.pivot
        const assemblySymmetryParams = AssemblySymmetryProvider.getParams(s.cell.obj?.data!);

        if (s.properties) {
            const oldParams = s.properties.cell.transform.params?.properties[AssemblySymmetryProvider.descriptor.name];
            if (PD.areEqual(assemblySymmetryParams, oldParams, values)) return;

            const b = this.plugin.state.data.build();
            b.to(s.properties.cell).update(old => {
                old.properties[AssemblySymmetryProvider.descriptor.name] = values;
            });
            await this.plugin.updateDataState(b);
        } else {
            const pd = this.plugin.customStructureProperties.getParams(s.cell.obj?.data);
            const params = PD.getDefaultValues(pd);
            if (PD.areEqual(assemblySymmetryParams, params.properties[AssemblySymmetryProvider.descriptor.name], values)) return;
            params.properties[AssemblySymmetryProvider.descriptor.name] = values;
            await this.plugin.builders.structure.insertStructureProperties(s.cell, params);
        }
    }

    paramsOnChange = (options: AssemblySymmetryProps) => {
        this.setState({ options })
        this.updateAssemblySymmetry(options)
    }

    get hasAssemblySymmetry3D() {
        return !!this.pivot.genericRepresentations?.filter(r => r.cell.transform.transformer.id === AssemblySymmetry3D.id)[0]
    }

    get enable() {
        return !this.hasAssemblySymmetry3D && this.state.options.symmetryIndex !== -1
    }

    get noSymmetries() {
        const structure = this.pivot.cell.obj?.data
        const data = structure && AssemblySymmetryDataProvider.get(structure).value
        return data && data.filter(sym => sym.symbol !== 'C1').length === 0
    }

    renderParams() {
        return <>
            <ParameterControls params={this.params} values={this.state.options} onChangeValues={this.paramsOnChange} />
        </>;
    }

    renderControls() {
        if (!this.pivot) return null;
        if (this.noSymmetries) return this.renderNoSymmetries();
        if (this.enable) return this.renderEnable();
        return this.renderParams();
    }
}