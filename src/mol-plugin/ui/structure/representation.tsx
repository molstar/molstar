/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { PluginUIComponent } from '../base';
import { Structure, StructureElement } from '../../../mol-model/structure';
import { isEmptyLoci } from '../../../mol-model/loci';
import { ColorOptions } from '../controls/parameters';
import { Color } from '../../../mol-util/color';
import { ButtonSelect, Options } from '../controls/common';
import { StructureSelectionQueries as Q } from '../../util/structure-selection-helper';

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
        this.plugin.helpers.structureRepresentation.set('remove', value, this.lociGetter)
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
                        <option key={-1} value={-1}>TODO: All</option>
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
        const { structureSelection: sel, structureRepresentation: rep } = this.plugin.helpers
        const lociGetter = (structure: Structure) => {
            const loci = this.plugin.helpers.structureSelectionManager.get(structure)
            return isEmptyLoci(loci) ? StructureElement.Loci.none(structure) : loci
        }

        sel.set('add', Q.all())
        await rep.set('add', 'cartoon', lociGetter)
        await rep.set('add', 'carbohydrate', lociGetter)

        sel.set('only', Q.ligandsPlusConnected())
        sel.set('add', Q.branchedConnectedOnly())
        sel.set('add', Q.water())
        await rep.set('add', 'ball-and-stick', lociGetter)

        sel.set('remove', Q.all())
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
        </div>
    }
}