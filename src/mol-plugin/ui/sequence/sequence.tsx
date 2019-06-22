/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react'
import { PluginUIComponent } from '../base';
import { Interactivity } from '../../util/interactivity';
import { MarkerAction } from '../../../mol-util/marker-action';
import { ButtonsType, ModifiersKeys, getButtons, getModifiers } from '../../../mol-util/input/input-observer';
import { ValueBox } from '../../../mol-util';
import { Residue } from './residue';
import { SequenceWrapper } from './wrapper';

type SequenceProps = { sequenceWrapper: SequenceWrapper.Any }
type SequenceState = { markerData: ValueBox<Uint8Array> }

function getState(markerData: ValueBox<Uint8Array>) {
    return { markerData: ValueBox.withValue(markerData, markerData.value) }
}

// TODO: this is really inefficient and should be done using a canvas.
export class Sequence<P extends SequenceProps> extends PluginUIComponent<P, SequenceState> {
    state = {
        markerData: ValueBox.create(this.props.sequenceWrapper.markerArray)
    }

    private setMarkerData(markerData: ValueBox<Uint8Array>) {
        this.setState(getState(markerData))
    }

    private lociHighlightProvider = (loci: Interactivity.Loci, action: MarkerAction) => {
        const changed = this.props.sequenceWrapper.markResidue(loci.loci, action)
        if (changed) this.setMarkerData(this.state.markerData)
    }

    private lociSelectionProvider = (loci: Interactivity.Loci, action: MarkerAction) => {
        const changed = this.props.sequenceWrapper.markResidue(loci.loci, action)
        if (changed) this.setMarkerData(this.state.markerData)
    }

    static getDerivedStateFromProps(nextProps: SequenceProps, prevState: SequenceState): SequenceState | null {
        if (prevState.markerData.value !== nextProps.sequenceWrapper.markerArray) {
            return getState(ValueBox.create(nextProps.sequenceWrapper.markerArray))
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

    highlight(seqId?: number, modifiers?: ModifiersKeys) {
        const ev = { current: Interactivity.Loci.Empty, modifiers }
        if (seqId !== undefined) {
            const loci = this.props.sequenceWrapper.getLoci(seqId);
            if (loci.elements.length > 0) ev.current = { loci };
        }
        this.plugin.behaviors.interaction.highlight.next(ev)
    }

    click(seqId: number | undefined, buttons: ButtonsType, modifiers: ModifiersKeys) {
        const ev = { current: Interactivity.Loci.Empty, buttons, modifiers }
        if (seqId !== undefined) {
            const loci = this.props.sequenceWrapper.getLoci(seqId);
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
        const sw = this.props.sequenceWrapper

        const elems: JSX.Element[] = [];
        for (let i = 0, il = sw.length; i < il; ++i) {
            elems[elems.length] = <Residue
                seqIdx={i}
                label={sw.residueLabel(i)}
                parent={this}
                marker={markerData.value[i]}
                color={sw.residueColor(i)}
                key={i}
            />;
        }

        return <div
            className='msp-sequence-wrapper'
            onContextMenu={this.contextMenu}
            onMouseDown={this.mouseDown}
        >
            {elems}
        </div>;
    }
}
