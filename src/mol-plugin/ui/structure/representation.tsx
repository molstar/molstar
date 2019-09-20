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
import { ButtonSelect, Options } from '../controls/common'
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualQuality, VisualQualityOptions } from '../../../mol-geo/geometry/base';
import { StructureRepresentationPresets as P } from '../../util/structure-representation-helper';
import { camelCaseToWords } from '../../../mol-util/string';
import { CollapsableControls } from '../base';

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
                    <optgroup label='Show'>{Options(types)}</optgroup>
                </ButtonSelect>
                <ButtonSelect label='Hide' onChange={this.hide}>
                    <optgroup label='Clear'>
                        <option key={'__all__'} value={'__all__'}>All</option>
                    </optgroup>
                    <optgroup label='Hide'>{Options(types)}</optgroup>
                </ButtonSelect>
                <ButtonSelect label='Color' onChange={this.color}>
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
        const loci = this.plugin.helpers.structureSelectionManager.get(structure)
        return isEmptyLoci(loci) ? StructureElement.Loci.none(structure) : loci
    }
}

export class StructureRepresentationControls extends CollapsableControls {
    preset = async (value: string) => {
        const presetFn = P[value as keyof typeof P]
        if (presetFn) {
            await presetFn(this.plugin.helpers.structureRepresentation)
        }
    }

    onChange = async (p: { param: PD.Base<any>, name: string, value: any }) => {
        if (p.name === 'options') {
            await this.plugin.helpers.structureRepresentation.setIgnoreHydrogens(!p.value.showHydrogens)
            await this.plugin.helpers.structureRepresentation.setQuality(p.value.visualQuality)
            this.forceUpdate()
        }
    }

    get params () {
        const { options } = this.values
        return {
            options: PD.Group({
                showHydrogens: PD.Boolean(options.showHydrogens),
                visualQuality: PD.Select<VisualQuality>(options.visualQuality, VisualQualityOptions),
            }, { isExpanded: true })
        }
    }

    get values () {
        const { structureRepresentation: rep } = this.plugin.helpers
        return {
            options: {
                showHydrogens: !rep.ignoreHydrogens,
                visualQuality: rep.quality,
            }
        }
    }

    defaultState() {
        return {
            isCollapsed: false,
            header: 'Representation'
        }
    }

    renderControls() {
        const presets = Object.keys(P).map(name => {
            return [name, camelCaseToWords(name)] as [string, string]
        })

        return <div>
            <div className='msp-control-row'>
                <div className='msp-select-row'>
                    <ButtonSelect label='Preset' onChange={this.preset}>
                        <optgroup label='Preset'>{Options(presets)}</optgroup>
                    </ButtonSelect>
                </div>
            </div>
            <EverythingStructureRepresentationControls />
            <SelectionStructureRepresentationControls />

            <ParameterControls params={this.params} values={this.values} onChange={this.onChange} />
        </div>
    }
}