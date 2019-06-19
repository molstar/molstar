/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react'
import { StructureSelection, StructureQuery } from '../../../mol-model/structure';
import { PluginUIComponent } from '../base';
import { Interactivity } from '../../util/interactivity';
import { MarkerAction } from '../../../mol-util/marker-action';
import { ButtonsType, ModifiersKeys, getButtons, getModifiers } from '../../../mol-util/input/input-observer';
import { ValueBox } from '../../../mol-util';
import { createResidueQuery, markResidue, StructureSeq } from './util';
import { Residue } from './residue';

type BaseSequenceProps = { structureSeq: StructureSeq, markerArray: Uint8Array }
type BaseSequenceState = { markerData: ValueBox<Uint8Array> }

// TODO: this is really inefficient and should be done using a canvas.
export class BaseSequence extends PluginUIComponent<BaseSequenceProps, BaseSequenceState> {
    state = {
        markerData: ValueBox.create(new Uint8Array(this.props.markerArray))
    }

    private lociHighlightProvider = (loci: Interactivity.Loci, action: MarkerAction) => {
        const { markerData } = this.state;
        const changed = markResidue(loci.loci, this.props.structureSeq, markerData.value, action)
        if (changed) this.setState({ markerData: ValueBox.withValue(markerData, markerData.value) })
    }

    private lociSelectionProvider = (loci: Interactivity.Loci, action: MarkerAction) => {
        const { markerData } = this.state;
        const changed = markResidue(loci.loci, this.props.structureSeq, markerData.value, action)
        if (changed) this.setState({ markerData: ValueBox.withValue(markerData, markerData.value) })
    }

    static getDerivedStateFromProps(nextProps: BaseSequenceProps, prevState: BaseSequenceState): BaseSequenceState | null {
        if (prevState.markerData.value !== nextProps.markerArray) {
            return { markerData: ValueBox.create(nextProps.markerArray) }
        }
        return null
    }

    componentDidMount() {
        this.plugin.interactivity.lociHighlights.addProvider(this.lociHighlightProvider)
        this.plugin.interactivity.lociSelections.addProvider(this.lociSelectionProvider)
    }

    componentWillUnmount() {
        this.plugin.interactivity.lociHighlights.removeProvider(this.lociHighlightProvider)
        this.plugin.interactivity.lociSelections.removeProvider(this.lociSelectionProvider)
    }

    getLoci(seqId: number) {
        const { structure, seq } = this.props.structureSeq
        const query = createResidueQuery(seq.entityId, seqId);
        return StructureSelection.toLoci2(StructureQuery.run(query, structure));
    }

    highlight(seqId?: number, modifiers?: ModifiersKeys) {
        const ev = { current: Interactivity.Loci.Empty, modifiers }
        if (seqId !== undefined) {
            const loci = this.getLoci(seqId);
            if (loci.elements.length > 0) ev.current = { loci };
        }
        this.plugin.behaviors.interaction.highlight.next(ev)
    }

    click(seqId: number | undefined, buttons: ButtonsType, modifiers: ModifiersKeys) {
        const ev = { current: Interactivity.Loci.Empty, buttons, modifiers }
        if (seqId !== undefined) {
            const loci = this.getLoci(seqId);
            if (loci.elements.length > 0) ev.current = { loci };
        }
        this.plugin.behaviors.interaction.click.next(ev)
    }

    contextMenu = (e: React.MouseEvent) => {
        e.preventDefault()
    }

    mouseDown = (e: React.MouseEvent) => {
        const buttons = getButtons(e.nativeEvent)
        const modifiers = getModifiers(e.nativeEvent)
        this.click(undefined, buttons, modifiers);
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
