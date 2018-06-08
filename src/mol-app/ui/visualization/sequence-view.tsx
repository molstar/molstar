/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react'
import { View } from '../view';
import { SequenceViewController } from '../../controller/visualization/sequence-view';
import { Structure, StructureSequence, Queries, Selection } from 'mol-model/structure';
import { Context } from '../../context/context';
import { InteractivityEvents } from '../../event/basic';
import { SyncRuntimeContext } from 'mol-task/execution/synchronous';

export class SequenceView extends View<SequenceViewController, {}, {}> {
    render() {
        const s = this.controller.latestState.structure;
        if (!s) return <div className='molstar-sequence-view-wrap'>No structure available.</div>;

        const seqs = Structure.getModels(s)[0].sequence.sequences;
        return <div className='molstar-sequence-view-wrap'>
            {seqs.map((seq, i) => <EntitySequence key={i} ctx={this.controller.context} seq={seq} structure={s} /> )}
        </div>;
    }
}

function createQuery(entityId: string, label_seq_id: number) {
    return Queries.generators.atoms({
        entityTest: l => Queries.props.entity.id(l) === entityId,
        residueTest: l => Queries.props.residue.label_seq_id(l) === label_seq_id
    });
}

// TODO: this is really ineffective and should be done using a canvas.
class EntitySequence extends React.Component<{ ctx: Context, seq: StructureSequence.Entity, structure: Structure }> {

    async raiseInteractityEvent(seqId?: number) {
        if (typeof seqId === 'undefined') {
            InteractivityEvents.HighlightElementLoci.dispatch(this.props.ctx, void 0);
            return;
        }

        const query = createQuery(this.props.seq.entityId, seqId);
        const loci = Selection.toLoci(await query(this.props.structure, SyncRuntimeContext));
        InteractivityEvents.HighlightElementLoci.dispatch(this.props.ctx, loci);
    }


    render() {
        const { ctx, seq } = this.props;
        const { offset, sequence } = seq.sequence;

        const elems: JSX.Element[] = [];
        for (let i = 0, _i = sequence.length; i < _i; i++) {
            elems[elems.length] = <ResidueView ctx={ctx} seqId={offset + i} letter={sequence[i]} parent={this} key={i} />;
        }

        return <div style={{ wordWrap: 'break-word' }}>
            <span style={{ fontWeight: 'bold' }}>{this.props.seq.entityId}:{offset}&nbsp;</span>
            {elems}
        </div>;
    }
}

class ResidueView extends React.Component<{ ctx: Context, seqId: number, letter: string, parent: EntitySequence }, { isHighlighted: boolean }> {
    state = { isHighlighted: false }

    mouseEnter = () => {
        this.setState({ isHighlighted: true });
        this.props.parent.raiseInteractityEvent(this.props.seqId);
    }

    mouseLeave = () => {
        this.setState({ isHighlighted: false });
        this.props.parent.raiseInteractityEvent();
    }

    render() {
        return <span onMouseEnter={this.mouseEnter} onMouseLeave={this.mouseLeave}
            style={{ cursor: 'pointer', backgroundColor: this.state.isHighlighted ? 'yellow' : void 0 }}>
            {this.props.letter}
        </span>;
    }
}