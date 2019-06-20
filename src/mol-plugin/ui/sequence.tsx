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
import { Structure } from '../../mol-model/structure';
import { SequenceWrapper } from './sequence/util';
import { PolymerSequenceWrapper } from './sequence/polymer';
import { StructureElementSelectionManager } from '../util/structure-element-selection';
import { MarkerAction } from '../../mol-util/marker-action';

function getSequenceWrappersForStructure(structure: Structure, structureSelection: StructureElementSelectionManager) {
    const sequenceWrappers: SequenceWrapper.Any[] = []

    structure.units.forEach(unit => {
        if (unit.polymerElements.length === 0) return

        const sw = new PolymerSequenceWrapper({ structure, unit })
        sw.markResidue(structureSelection.get(structure), MarkerAction.Select)
        sequenceWrappers.push(sw)
    })

    return sequenceWrappers
}

export class SequenceView extends PluginUIComponent<{ }, { }> {
    private spine: StateTreeSpine.Impl

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

    render() {
        const structure = this.getStructure();
        if (!structure) return <div className='msp-sequence'>
            <div className='msp-sequence-entity'>No structure available</div>
        </div>;

        const { structureSelection } = this.plugin.helpers
        const sequenceWrappers = getSequenceWrappersForStructure(structure, structureSelection)
        return <div className='msp-sequence'>
            {sequenceWrappers.map((sequenceWrapper, i) => {
                return <Sequence key={i} sequenceWrapper={sequenceWrapper} />
            })}
        </div>;
    }
}