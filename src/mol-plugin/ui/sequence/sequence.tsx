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
import { SequenceWrapper } from './wrapper';
import { StructureElement, StructureProperties, Unit } from '../../../mol-model/structure';
import { Subject } from 'rxjs';
import { debounceTime } from 'rxjs/operators';
import { Color } from '../../../mol-util/color';

type SequenceProps = {
    sequenceWrapper: SequenceWrapper.Any,
    sequenceNumberPeriod?: number,
    hideSequenceNumbers?: boolean
}

/** Note, if this is changed, the CSS for `msp-sequence-number` needs adjustment too */
const MaxSequenceNumberSize = 5

// TODO: this is somewhat inefficient and should be done using a canvas.
export class Sequence<P extends SequenceProps> extends PluginUIComponent<P> {
    private parentDiv = React.createRef<HTMLDivElement>();
    private lastMouseOverSeqIdx = -1;
    private highlightQueue = new Subject<{ seqIdx: number, buttons: number, modifiers: ModifiersKeys }>();

    private lociHighlightProvider = (loci: Interactivity.Loci, action: MarkerAction) => {
        const changed = this.props.sequenceWrapper.markResidue(loci.loci, action)
        if (changed) this.updateMarker();
    }

    private lociSelectionProvider = (loci: Interactivity.Loci, action: MarkerAction) => {
        const changed = this.props.sequenceWrapper.markResidue(loci.loci, action)
        if (changed) this.updateMarker();
    }

    private get sequenceNumberPeriod() {
        if (this.props.sequenceNumberPeriod !== undefined) {
            return this.props.sequenceNumberPeriod as number
        }
        if (this.props.sequenceWrapper.length > 10) return 10
        const lastSeqNum = this.getSequenceNumber(this.props.sequenceWrapper.length - 1)
        if (lastSeqNum.length > 1) return 5
        return 1
    }

    componentDidMount() {
        this.plugin.interactivity.lociHighlights.addProvider(this.lociHighlightProvider)
        this.plugin.interactivity.lociSelects.addProvider(this.lociSelectionProvider)

        this.subscribe(debounceTime<{ seqIdx: number, buttons: number, modifiers: ModifiersKeys }>(15)(this.highlightQueue), (e) => {
            this.hover(e.seqIdx < 0 ? void 0 : e.seqIdx, e.buttons, e.modifiers);
        });

        // this.updateMarker()
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

    private getResidueClass(seqIdx: number, label: string) {
        return label.length > 1
            ? (seqIdx === 0 ? 'msp-sequence-residue-long-begin' : 'msp-sequence-residue-long')
            : void 0
    }

    private residue(seqIdx: number, label: string, marker: number, color: Color) {
        return <span key={seqIdx} data-seqid={seqIdx} style={{ color: Color.toStyle(color), backgroundColor: this.getBackgroundColor(marker) }} className={this.getResidueClass(seqIdx, label)}>{label}</span>;
    }

    private getSequenceNumberClass(seqIdx: number, seqNum: string, label: string) {
        const classList = ['msp-sequence-number']
        if (seqNum.startsWith('-')) {
            if (label.length > 1 && seqIdx > 0) classList.push('msp-sequence-number-long-negative')
            else classList.push('msp-sequence-number-negative')
        } else {
            if (label.length > 1 && seqIdx > 0) classList.push('msp-sequence-number-long')
        }
        return classList.join(' ')
    }

    private location = StructureElement.Location.create();
    private getSequenceNumber(seqIdx: number) {
        let seqNum = ''
        const loci = this.props.sequenceWrapper.getLoci(seqIdx)
        const l = StructureElement.Loci.getFirstLocation(loci, this.location);
        if (l) {
            if (Unit.isAtomic(l.unit)) {
                const seqId = StructureProperties.residue.auth_seq_id(l)
                const insCode = StructureProperties.residue.pdbx_PDB_ins_code(l)
                seqNum = `${seqId}${insCode ? insCode : ''}`
            } else if (Unit.isCoarse(l.unit)) {
                seqNum = `${seqIdx + 1}`
            }
        }
        return seqNum
    }

    private padSeqNum(n: string) {
        if (n.length < MaxSequenceNumberSize) return n + new Array(MaxSequenceNumberSize - n.length + 1).join('\u00A0');
        return n;
    }
    private getSequenceNumberSpan(seqIdx: number, label: string) {
        const seqNum = this.getSequenceNumber(seqIdx)
        return <span key={`marker-${seqIdx}`} className={this.getSequenceNumberClass(seqIdx, seqNum, label)}>{this.padSeqNum(seqNum)}</span>
    }

    private updateMarker() {
        if (!this.parentDiv.current) return;
        const xs = this.parentDiv.current.children;
        const { markerArray } = this.props.sequenceWrapper;
        const hasNumbers = !this.props.hideSequenceNumbers, period = this.sequenceNumberPeriod;

        // let first: HTMLSpanElement | undefined;

        let o = 0;
        for (let i = 0, il = markerArray.length; i < il; i++) {
            if (hasNumbers && i % period === 0 && i < il) o++;
            // o + 1 to account for help icon
            const span = xs[o] as HTMLSpanElement;
            if (!span) return;
            o++;

            // if (!first && markerArray[i] > 0) {
            //     first = span;
            // }

            const backgroundColor = this.getBackgroundColor(markerArray[i]);
            if (span.style.backgroundColor !== backgroundColor) span.style.backgroundColor = backgroundColor;
        }

        // if (first) {
        //     first.scrollIntoView({ block: 'nearest' });
        // }
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
        const sw = this.props.sequenceWrapper

        const elems: JSX.Element[] = [];

        const hasNumbers = !this.props.hideSequenceNumbers, period = this.sequenceNumberPeriod;
        for (let i = 0, il = sw.length; i < il; ++i) {
            const label = sw.residueLabel(i)
            // add sequence number before name so the html element do not get separated by a line-break
            if (hasNumbers && i % period === 0 && i < il) {
                elems[elems.length] = this.getSequenceNumberSpan(i, label)
            }
            elems[elems.length] = this.residue(i, label, sw.markerArray[i], sw.residueColor(i));
        }

        // calling .updateMarker here is neccesary to ensure existing
        // residue spans are updated as react won't update them
        this.updateMarker()

        // <span className={`msp-icon msp-icon-help`} style={{ cursor: 'help' }} title='This shows a single sequence. Use the menu on the left to show a different sequence.' />

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
