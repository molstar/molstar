/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { PluginUIComponent } from '../base';
import { Structure, StructureElement } from '../../../mol-model/structure';
import { isEmptyLoci } from '../../../mol-model/loci';
import { ColorOptions, ParameterControls } from '../controls/parameters';
import { Color } from '../../../mol-util/color';
import { ButtonSelect, Options } from '../controls/common';
import { StructureSelectionQueries as Q } from '../../util/structure-selection-helper';
import { MolScriptBuilder as MS } from '../../../mol-script/language/builder';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';

abstract class BaseStructureRepresentationControls extends PluginUIComponent {
    onChange = (value: string) => {
        console.log('onChange', value)
    }

    abstract label: string
    abstract lociGetter(structure: Structure): StructureElement.Loci

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
        const { types } = this.plugin.structureRepresentation.registry

        return <div className='msp-control-row'>
            <span title={this.label}>{this.label}</span>
            <div className='msp-select-row'>
                <ButtonSelect label='Show' onChange={this.show}>
                    <optgroup label='Show'>
                        {Options(types)}
                    </optgroup>
                </ButtonSelect>
                <ButtonSelect label='Hide' onChange={this.hide}>
                    <optgroup label='Clear'>
                        <option key={'__all__'} value={'__all__'}>All</option>
                    </optgroup>
                    <optgroup label='Hide'>
                        {Options(types)}
                    </optgroup>
                </ButtonSelect>
                <ButtonSelect label='Color' onChange={this.color}>
                    <optgroup label='Clear'>
                        <option key={-1} value={-1}>Theme</option>
                    </optgroup>
                    <optgroup label='Color'>
                        {ColorOptions()}
                    </optgroup>
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
        const loci = this.plugin.helpers.structureSelectionManager.get(structure)
        return isEmptyLoci(loci) ? StructureElement.Loci.none(structure) : loci
    }
}

export class StructureRepresentationControls extends PluginUIComponent {
    preset = async () => {
        const { structureRepresentation: rep } = this.plugin.helpers
        await rep.clear()
        await rep.setFromExpression('add', 'cartoon', Q.all)
        await rep.setFromExpression('add', 'carbohydrate', Q.all)
        await rep.setFromExpression('add', 'ball-and-stick', MS.struct.modifier.union([
            MS.struct.combinator.merge([ Q.ligandsPlusConnected, Q.branchedConnectedOnly, Q.water ])
        ]))
    }

    onChange = async (p: { param: PD.Base<any>, name: string, value: any }) => {
        if (p.name === 'showHydrogens') {
            await this.plugin.helpers.structureRepresentation.setIgnoreHydrogens(!p.value)
        }
        this.forceUpdate()
    }

    get params () {
        const values = this.values
        return {
            showHydrogens: PD.Boolean(values.showHydrogens)
        }
    }

    get values () {
        const { structureRepresentation: rep } = this.plugin.helpers
        return {
            showHydrogens: !rep.ignoreHydrogens
        }
    }

    render() {
        return <div className='msp-transform-wrapper'>
            <div className='msp-transform-header'>
                <button className='msp-btn msp-btn-block'>Representation</button>
            </div>
            <div className='msp-btn-row-group'>
                <button className='msp-btn msp-btn-block msp-form-control' onClick={() => this.preset()}>Preset</button>
            </div>
            <EverythingStructureRepresentationControls />
            <SelectionStructureRepresentationControls />

            <ParameterControls params={this.params} values={this.values} onChange={this.onChange} />
        </div>
    }
}