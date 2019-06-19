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
import { MarkerAction } from '../../mol-util/marker-action';
import { PolymerSequence } from './sequence/polymer';
import { StructureSeq, markResidue } from './sequence/util';

function getStructureSeqKey(structureSeq: StructureSeq) {
    const { structure, seq } = structureSeq
    const strucHash = structure.parent ? structure.parent.hashCode : structure.hashCode
    return `${strucHash}|${seq.entityId}`
}

export class SequenceView extends PluginUIComponent<{ }, { }> {
    private spine: StateTreeSpine.Impl
    private markerArrays = new Map<string, Uint8Array>()

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

    private getMarkerArray(structureSeq: StructureSeq): Uint8Array {
        const { structure, seq } = structureSeq
        const key = getStructureSeqKey(structureSeq)
        let markerArray = this.markerArrays.get(key)
        if (!markerArray) {
            markerArray = new Uint8Array(seq.sequence.sequence.length)
            this.markerArrays.set(key, markerArray)
        }
        const loci = this.plugin.helpers.structureSelection.get(structure)
        markerArray.fill(0)
        markResidue(loci, structureSeq, markerArray, MarkerAction.Select)
        return markerArray
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

        const seqs = structure.models[0].sequence.sequences;
        return <div className='msp-sequence'>
            {seqs.map((seq, i) => {
                const structureSeq = { structure, seq }
                const markerArray = this.getMarkerArray(structureSeq)
                return <PolymerSequence key={i} structureSeq={structureSeq} markerArray={markerArray} />
            })}
        </div>;
    }
}