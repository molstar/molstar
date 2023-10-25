/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Yakov Pechersky <ffxen158@gmail.com>
 */

import { Structure, StructureElement, Unit } from '../../mol-model/structure';
import { SecondaryStructureProvider } from '../../mol-model-props/computed/secondary-structure';
import { ModelSecondaryStructure } from '../../mol-model-formats/structure/property/secondary-structure';
import { Sequence } from './sequence';
import { SequenceWrapper } from './wrapper';

type SecondaryStructureSequenceProps = {
    sequenceWrapper: SequenceWrapper.Any,
    sequenceNumberPeriod?: number,
    hideSequenceNumbers?: boolean,
    hideSecondaryStructure?: boolean,
}

export class SecondaryStructureSequence extends Sequence<SecondaryStructureSequenceProps> {

    protected updateMarker() {
        if (!this.parentDiv.current) return;
        const xs = this.parentDiv.current.querySelectorAll('.msp-sequence-missing, .msp-sequence-present');
        const { markerArray } = this.props.sequenceWrapper;
        const secondarySpanDiv = this.parentDiv.current.querySelector('.msp-sequence-secondary');
        const emptySS = (!!secondarySpanDiv) && /^[\s\u200b]+$/.test(secondarySpanDiv.textContent ?? '');
        const overlays = [' '];

        for (let i = 0, il = markerArray.length; i < il; i++) {
            const span = xs[i] as HTMLSpanElement | undefined;
            if (!span) continue;

            const backgroundColor = this.getBackgroundColor(markerArray[i]);
            if (span.style.backgroundColor !== backgroundColor) span.style.backgroundColor = backgroundColor;
            if (emptySS) {
                // only calculate if needed
                const ss = this.secondaryStructureKind(i);
                overlays.push(ss[0]);
            }
        }
        if (emptySS) {
            secondarySpanDiv.textContent = overlays.join('\u200b');
        }

    }

    protected getSequenceNumberClass(seqIdx: number, seqNum: string, label: string) {
        const suffix = this.props.hideSecondaryStructure ? '' : ' msp-sequence-number-above';
        return super.getSequenceNumberClass(seqIdx, seqNum, label) + suffix;
    }

    protected secondaryStructureKind(i: number): string {
        const loci = this.props.sequenceWrapper.getLoci(i);
        const l = StructureElement.Loci.getFirstLocation(loci, this.location);
        if (!l || !Unit.isAtomic(l.unit)) return ' ';
        const secStruc = SecondaryStructureProvider.get(l.structure).value?.get(l.unit.invariantId);
        if (!secStruc) return ' ';
        const elem =
        secStruc.elements[
            secStruc.key[secStruc.getIndex(l.unit.residueIndex[l.element])]
        ];
        if (!elem) return ' ';
        return elem.kind !== 'none' ? elem.kind : ' ';
    }


    render() {
        const sw = this.props.sequenceWrapper;
        const structure: Structure = this.props.sequenceWrapper.data.structure;
        const secondaryStructure = ModelSecondaryStructure.Provider.get(structure.model);

        const elems: JSX.Element[] = [];

        const hasNumbers = !this.props.hideSequenceNumbers, period = this.sequenceNumberPeriod;
        const overlays = [' '];
        for (let i = 0, il = sw.length; i < il; ++i) {
            const label = sw.residueLabel(i);
            // add sequence number before name so the html element do not get separated by a line-break
            if (hasNumbers && i % period === 0 && i < il) {
                elems[elems.length] = this.getSequenceNumberSpan(i, label);
            }
            elems[elems.length] = this.residue(i, label, sw.markerArray[i]);
            if (!this.props.hideSecondaryStructure) {
                if (!secondaryStructure) continue;
                const ss = this.secondaryStructureKind(i);
                overlays.push(ss[0]);
            }
        }
        elems[elems.length] = (
            <div className="msp-sequence-secondary">{overlays.join('\u200b')}</div>
        );

        // calling .updateMarker here is neccesary to ensure existing
        // residue spans are updated as react won't update them
        this.updateMarker();

        const className = this.props.hideSecondaryStructure ? 'msp-sequence-wrapper' : 'msp-sequence-wrapper msp-sequence-wrapper-secondary';

        return <div
            className={className}
            onContextMenu={this.contextMenu}
            onMouseDown={this.mouseDown}
            onMouseUp={this.mouseUp}
            onMouseMove={this.mouseMove}
            onMouseLeave={this.mouseLeave}
            ref={this.parentDiv}
        >
            {elems}
        </div>;
    }
}
