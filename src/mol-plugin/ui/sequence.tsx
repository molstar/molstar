/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react'
import { PluginUIComponent } from './base';
import { StateTreeSpine } from '../../mol-state/tree/spine';
import { PluginStateObject as SO } from '../state/objects';
import { Sequence } from './sequence/sequence';
import { Structure, StructureElement, StructureProperties as SP } from '../../mol-model/structure';
import { SequenceWrapper } from './sequence/util';
import { PolymerSequenceWrapper } from './sequence/polymer';
import { StructureElementSelectionManager } from '../util/structure-element-selection';
import { MarkerAction } from '../../mol-util/marker-action';
import { ParameterControls } from './controls/parameters';
import { ParamDefinition as PD } from '../../mol-util/param-definition';

function getSequenceWrapperForStructure(index: number, structure: Structure, structureSelection: StructureElementSelectionManager): SequenceWrapper.Any | undefined {
    let j = 0
    for (let i = 0, il = structure.units.length; i < il; ++i) {
        const unit = structure.units[i]
        if (unit.polymerElements.length === 0) continue
        if (j === index) {
            const sw = new PolymerSequenceWrapper({ structure, unit })
            sw.markResidue(structureSelection.get(structure), MarkerAction.Select)
            return sw
        }
        j += 1
    }
}

function getPolymerOptionsForStructure(structure: Structure) {
    const options: [number, string][] = []

    let i = 0
    structure.units.forEach(unit => {
        if (unit.polymerElements.length === 0) return

        const l = StructureElement.create(unit, unit.elements[0])
        const entityDescription = SP.entity.pdbx_description(l)
        const label_asym_id = SP.chain.label_asym_id(l)
        const label = `${label_asym_id}: ${entityDescription}`

        options.push([ i, label ])
        i += 1
    })

    return options
}

export class SequenceView extends PluginUIComponent<{ }, { polymer: number }> {
    private spine: StateTreeSpine.Impl

    state = { polymer: 0 }

    componentDidMount() {
        this.spine = new StateTreeSpine.Impl(this.plugin.state.dataState.cells);

        this.subscribe(this.plugin.state.behavior.currentObject, o => {
            const current = this.plugin.state.dataState.cells.get(o.ref)!;
            this.spine.current = current
            this.forceUpdate();
        });

        this.subscribe(this.plugin.events.state.object.updated, ({ ref, state }) => {
            const current = this.spine.current;
            if (!current || current.sourceRef !== ref || current.state !== state) return;
            this.forceUpdate();
        });
    }

    private getStructure() {
        const so = this.spine && this.spine.getRootOfType(SO.Molecule.Structure)
        return so && so.data
    }

    private getParams(structure: Structure) {
        return {
            polymer: PD.Select(0, getPolymerOptionsForStructure(structure))
        }
    }

    private setParamProps = (p: { param: PD.Base<any>, name: string, value: any }) => {
        if (p.name === 'polymer') {
            this.setState({ polymer: p.value })
        }
    }

    render() {
        const structure = this.getStructure();
        if (!structure) return <div className='msp-sequence'>
            <div className='msp-sequence-entity'>No structure available</div>
        </div>;

        const { structureSelection } = this.plugin.helpers
        const params = this.getParams(structure)
        const sequenceWrapper = getSequenceWrapperForStructure(this.state.polymer, structure, structureSelection)

        return <div className='msp-sequence'>
            <ParameterControls params={params} values={this.state} onChange={this.setParamProps} />
            {sequenceWrapper !== undefined
                ? <Sequence sequenceWrapper={sequenceWrapper} />
                : <div className='msp-sequence-entity'>No sequence available</div>}
        </div>;
    }
}