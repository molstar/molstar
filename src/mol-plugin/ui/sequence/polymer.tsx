/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react'
import { StructureSelection, StructureQuery } from '../../../mol-model/structure';
import { createResidueQuery } from './util';
import { Residue } from './residue';
import { BaseSequence } from './base';

export class PolymerSequence extends BaseSequence {
    getLoci(seqId: number) {
        const { structure, seq } = this.props.structureSeq
        const query = createResidueQuery(seq.entityId, seqId);
        return StructureSelection.toLoci2(StructureQuery.run(query, structure));
    }

    render() {
        const { markerData } = this.state;
        const { seq } = this.props.structureSeq;
        const { offset, sequence } = seq.sequence;

        const elems: JSX.Element[] = [];
        for (let i = 0, _i = sequence.length; i < _i; i++) {
            elems[elems.length] = <Residue seqId={offset + i + 1} letter={sequence[i]} parent={this} marker={markerData.value[i]} key={i} />;
        }

        return <div
            className='msp-sequence-entity'
            onContextMenu={this.contextMenu}
            onMouseDown={this.mouseDown}
        >
            <span style={{ fontWeight: 'bold' }}>{seq.entityId}:{offset}&nbsp;</span>
            {elems}
        </div>;
    }
}
