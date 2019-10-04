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
import { ButtonsType, ModifiersKeys, getButtons, getModifiers, MouseModifiers } from '../../../mol-util/input/input-observer';
import { SequenceWrapper } from './wrapper';
import { StructureElement } from '../../../mol-model/structure';
import { Subject } from 'rxjs';
import { debounceTime } from 'rxjs/operators';
import { Color } from '../../../mol-util/color';

type SequenceProps = { sequenceWrapper: SequenceWrapper.Any }

// TODO: this is somewhat inefficient and should be done using a canvas.
export class Sequence<P extends SequenceProps> extends PluginUIComponent<P> {
    private parentDiv = React.createRef<HTMLDivElement>();
    private lastMouseOverSeqIdx = -1;
    private highlightQueue = new Subject<{ seqIdx: number, buttons: number, modifiers: MouseModifiers }>();

    private lociHighlightProvider = (loci: Interactivity.Loci, action: MarkerAction) => {
        const changed = this.props.sequenceWrapper.markResidue(loci.loci, action)
        if (changed) this.updateMarker();
    }

    private lociSelectionProvider = (loci: Interactivity.Loci, action: MarkerAction) => {
        const changed = this.props.sequenceWrapper.markResidue(loci.loci, action)
        if (changed) this.updateMarker();
    }

    componentDidMount() {
        this.plugin.interactivity.lociHighlights.addProvider(this.lociHighlightProvider)
        this.plugin.interactivity.lociSelects.addProvider(this.lociSelectionProvider)

        this.subscribe(debounceTime<{ seqIdx: number, buttons: number, modifiers: MouseModifiers }>(15)(this.highlightQueue), (e) => {
            this.hover(e.seqIdx < 0 ? void 0 : e.seqIdx, e.buttons, e.modifiers);
        });
    }

    componentWillUnmount() {
        this.plugin.interactivity.lociHighlights.removeProvider(this.lociHighlightProvider)
        this.plugin.interactivity.lociSelects.removeProvider(this.lociSelectionProvider)
    }

    hover(seqId: number | undefined, buttons: ButtonsType, modifiers: ModifiersKeys) {
        const ev = { current: Interactivity.Loci.Empty, buttons, modifiers }
        if (seqId !== undefined) {
            const loci = this.props.sequenceWrapper.getLoci(seqId);
            if (!StructureElement.Loci.isEmpty(loci)) ev.current = { loci };
        }
        this.plugin.behaviors.interaction.hover.next(ev)
    }

    click(seqId: number | undefined, buttons: ButtonsType, modifiers: ModifiersKeys) {
        const ev = { current: Interactivity.Loci.Empty, buttons, modifiers }
        if (seqId !== undefined) {
            const loci = this.props.sequenceWrapper.getLoci(seqId);
            if (!StructureElement.Loci.isEmpty(loci)) ev.current = { loci };
        }
        this.plugin.behaviors.interaction.click.next(ev)
    }

    contextMenu = (e: React.MouseEvent) => {
        e.preventDefault()
    }

    mouseDown = (e: React.MouseEvent) => {
        e.stopPropagation();

        const buttons = getButtons(e.nativeEvent)
        const modifiers = getModifiers(e.nativeEvent)

        let seqIdx: number | undefined = undefined;
        const el = e.target as HTMLElement;
        if (el && el.getAttribute) {
            seqIdx = el.hasAttribute('data-seqid') ? +el.getAttribute('data-seqid')! : undefined;
        }
        this.click(seqIdx, buttons, modifiers);
    }

    private getBackgroundColor(marker: number) {
        // TODO: make marker color configurable
        if (typeof marker === 'undefined') console.error('unexpected marker value')
        return marker === 0 ? '' : marker % 2 === 0 ? 'rgb(51, 255, 25)' /* selected */ : 'rgb(255, 102, 153)' /* highlighted */;
    }

    private residue(seqIdx: number, label: string, marker: number, color: Color) {
        const margin = label.length > 1 ? (seqIdx === 0 ? `0px 2px 0px 0px` : `0px 2px 0px 2px`) : void 0
        return <span key={seqIdx} data-seqid={seqIdx} style={{ color: Color.toStyle(color), backgroundColor: this.getBackgroundColor(marker), margin }}>{label}</span>;
    }

    private updateMarker() {
        if (!this.parentDiv.current) return;
        const xs = this.parentDiv.current.children;
        const markerData = this.props.sequenceWrapper.markerArray;

        for (let i = 0, _i = markerData.length; i < _i; i++) {
            const span = xs[i] as HTMLSpanElement;
            if (!span) continue;

            const backgroundColor = this.getBackgroundColor(markerData[i]);
            if (span.style.backgroundColor !== backgroundColor) span.style.backgroundColor = backgroundColor;
        }
    }

    mouseMove = (e: React.MouseEvent) => {
        e.stopPropagation();

        const el = e.target as HTMLElement;
        if (!el || !el.getAttribute) {
            if (this.lastMouseOverSeqIdx === -1) return;

            this.lastMouseOverSeqIdx = -1;
            const buttons = getButtons(e.nativeEvent)
            const modifiers = getModifiers(e.nativeEvent)
            this.highlightQueue.next({ seqIdx: -1, buttons, modifiers })
            return;
        }
        const seqIdx = el.hasAttribute('data-seqid') ? +el.getAttribute('data-seqid')! : -1;
        if (this.lastMouseOverSeqIdx === seqIdx) {
            return;
        } else {
            const buttons = getButtons(e.nativeEvent)
            const modifiers = getModifiers(e.nativeEvent)
            this.lastMouseOverSeqIdx = seqIdx;
            this.highlightQueue.next({ seqIdx, buttons, modifiers })
        }
    }

    mouseLeave = (e: React.MouseEvent) => {
        if (this.lastMouseOverSeqIdx === -1) return;
        this.lastMouseOverSeqIdx = -1;
        const buttons = getButtons(e.nativeEvent)
        const modifiers = getModifiers(e.nativeEvent)
        this.highlightQueue.next({ seqIdx: -1, buttons, modifiers })
    }

    render() {
        const markerData = this.props.sequenceWrapper.markerArray;
        const sw = this.props.sequenceWrapper

        const elems: JSX.Element[] = [];
        for (let i = 0, il = sw.length; i < il; ++i) {
            elems[elems.length] = this.residue(i, sw.residueLabel(i), markerData[i], sw.residueColor(i));
            // TODO: add seq idx markers every N residues? Would need to modify "updateMarker"
        }

        return <div
            className='msp-sequence-wrapper msp-sequence-wrapper-non-empty'
            onContextMenu={this.contextMenu}
            onMouseDown={this.mouseDown}
            onMouseMove={this.mouseMove}
            onMouseLeave={this.mouseLeave}
            ref={this.parentDiv}
        >
            {elems}
        </div>;
    }
}
