/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { PluginUIComponent, CollapsableState, CollapsableProps } from '../base';
import { Structure, StructureElement } from '../../mol-model/structure';
import { isEmptyLoci } from '../../mol-model/loci';
import { ParameterControls } from '../controls/parameters';
import { Color } from '../../mol-util/color';
import { ButtonSelect, Options } from '../controls/common'
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { VisualQuality, VisualQualityOptions } from '../../mol-geo/geometry/base';
import { CollapsableControls } from '../base';
import { StateSelection, StateObject } from '../../mol-state';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { ColorOptions } from '../controls/color';
import { InteractionsProvider } from '../../mol-model-props/computed/interactions';

interface BaseStructureRepresentationControlsState {
    isDisabled: boolean
}

abstract class BaseStructureRepresentationControls extends PluginUIComponent<{}, BaseStructureRepresentationControlsState> {
    abstract label: string
    abstract lociGetter(structure: Structure): StructureElement.Loci

    state = {
        isDisabled: false
    }

    /** root structures */
    protected get structures() {
        return this.plugin.state.dataState.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Structure)).map(s => s.obj!.data)
    }

    /** applicable types */
    private get types() {
        const types: [string, string][] = []
        const structures = this.structures
        for (const e of this.plugin.structureRepresentation.registry.list) {
            if (structures.some(s => e.provider.isApplicable(s))) {
                types.push([e.name, e.provider.label])
            }
        }
        return types
    }

    private forceUpdateIfStructure = (obj?: StateObject) => {
        if (obj && obj.type === PluginStateObject.Molecule.Structure.type) {
            this.forceUpdate()
        }
    }

    componentDidMount() {
        this.subscribe(this.plugin.events.state.object.created, ({ obj }) => this.forceUpdateIfStructure(obj))

        this.subscribe(this.plugin.events.state.object.removed, ({ obj }) => this.forceUpdateIfStructure(obj))

        this.subscribe(this.plugin.state.dataState.events.isUpdating, v => this.setState({ isDisabled: v }))
    }

    show = (value: string) => {
        this.plugin.helpers.structureRepresentation.set('add', value, this.lociGetter)
    }

    hide = (value: string) => {
        if (value === '__all__') {
            const { types } = this.plugin.structureRepresentation.registry
            for (let i = 0, il = types.length; i < il; ++i) {
                this.plugin.helpers.structureRepresentation.set('remove', types[i][0], this.lociGetter)
            }
        } else {
            this.plugin.helpers.structureRepresentation.set('remove', value, this.lociGetter)
        }
    }

    color = (value: string) => {
        const color = Color(parseInt(value))
        this.plugin.helpers.structureOverpaint.set(color, this.lociGetter)
    }

    render() {
        const TypeOptions = Options(this.types)

        return <div className='msp-control-row'>
            <span title={this.label}>{this.label}</span>
            <div className='msp-select-row'>
                <ButtonSelect label='Show' onChange={this.show} disabled={this.state.isDisabled}>
                    <optgroup label='Show'>{TypeOptions}</optgroup>
                </ButtonSelect>
                <ButtonSelect label='Hide' onChange={this.hide} disabled={this.state.isDisabled}>
                    <optgroup label='Clear'>
                        <option key={'__all__'} value={'__all__'}>All</option>
                    </optgroup>
                    <optgroup label='Hide'>{TypeOptions}</optgroup>
                </ButtonSelect>
                <ButtonSelect label='Color' onChange={this.color} disabled={this.state.isDisabled}>
                    <optgroup label='Clear'>
                        <option key={-1} value={-1}>Theme</option>
                    </optgroup>
                    <optgroup label='Color'>{ColorOptions()}</optgroup>
                </ButtonSelect>
            </div>
        </div>
    }
}

class EverythingStructureRepresentationControls extends BaseStructureRepresentationControls {
    label = 'Everything'
    lociGetter = (structure: Structure) => {
        return StructureElement.Loci.all(structure)
    }
}

class SelectionStructureRepresentationControls extends BaseStructureRepresentationControls {
    label = 'Selection'
    lociGetter = (structure: Structure) => {
        const loci = this.plugin.managers.structure.selection.getLoci(structure)
        return isEmptyLoci(loci) ? StructureElement.Loci.none(structure) : loci
    }
}

interface StructureRepresentationControlsState extends CollapsableState {
    isDisabled: boolean
}

export class StructureRepresentationControls extends CollapsableControls<CollapsableProps, StructureRepresentationControlsState> {
    componentDidMount() {
        this.subscribe(this.plugin.state.dataState.events.isUpdating, v => this.setState({ isDisabled: v }))
    }

    onChange = async (p: { param: PD.Base<any>, name: string, value: any }) => {
        if (p.name === 'options') {
            await this.plugin.helpers.structureRepresentation.setIgnoreHydrogens(!p.value.showHydrogens)
            await this.plugin.helpers.structureRepresentation.setQuality(p.value.visualQuality)
            await this.plugin.helpers.structureRepresentation.setInteractionsProps(p.value.interactions)
            this.forceUpdate()
        }
    }

    get params () {
        const { options } = this.values
        return {
            options: PD.Group({
                showHydrogens: PD.Boolean(options.showHydrogens, { description: 'Toggle display of hydrogen atoms in representations' }),
                visualQuality: PD.Select<VisualQuality>(options.visualQuality, VisualQualityOptions, { description: 'Control the visual/rendering quality of representations' }),
                interactions: PD.Group(InteractionsProvider.defaultParams, { label: 'Non-covalent Interactions' }),
            }, { isExpanded: true })
        }
    }

    get values () {
        const { structureRepresentation: rep } = this.plugin.helpers
        return {
            options: {
                showHydrogens: !rep.ignoreHydrogens,
                visualQuality: rep.quality,
                interactions: rep.interactionProps,
            }
        }
    }

    defaultState() {
        return {
            isCollapsed: false,
            header: 'Representation',

            isDisabled: false
        }
    }

    renderControls() {
        return <div>
            <EverythingStructureRepresentationControls />
            <SelectionStructureRepresentationControls />

            <ParameterControls params={this.params} values={this.values} onChange={this.onChange} isDisabled={this.state.isDisabled} />
        </div>
    }
}