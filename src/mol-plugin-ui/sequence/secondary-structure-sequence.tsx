/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Yakov Pechersky <ffxen158@gmail.com>
 */

import { Structure, StructureElement, StructureProperties } from '../../mol-model/structure';
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
        const xs = this.parentDiv.current.querySelectorAll('.msp-sequence-present');
        const { markerArray } = this.props.sequenceWrapper;

        for (let i = 0, il = markerArray.length; i < il; i++) {
            const span = xs[i] as HTMLSpanElement | undefined;
            if (!span) continue;

            const backgroundColor = this.getBackgroundColor(markerArray[i]);
            if (span.style.backgroundColor !== backgroundColor) span.style.backgroundColor = backgroundColor;
        }

    }

    protected getSequenceNumberClass(seqIdx: number, seqNum: string, label: string) {
        const suffix = this.props.hideSecondaryStructure ? '' : ' msp-sequence-number-above';
        return super.getSequenceNumberClass(seqIdx, seqNum, label) + suffix;
    }

    protected getSequenceSecondaryStructureSpan(seqIdx: number) {
        const loci = this.props.sequenceWrapper.getLoci(seqIdx);
        const location = StructureElement.Loci.getFirstLocation(loci, this.location);
        if (!location) return;
        const structure: Structure = this.props.sequenceWrapper.data.structure;
        const secondaryStructure = ModelSecondaryStructure.Provider.get(structure.model);
        if (!secondaryStructure) return;
        const { kind } =
        secondaryStructure.elements[
            secondaryStructure.key[
                secondaryStructure.getIndex(
                    StructureProperties.residue.key(location)
                )
            ]
        ];
        if (kind !== 'none') return;
        const span = <span className="msp-sequence-secondary">{`\u200b${kind[0]}\u200b`}</span>;
        return span;
    }


    render() {
        const sw = this.props.sequenceWrapper;

        const elems: JSX.Element[] = [];

        const hasNumbers = !this.props.hideSequenceNumbers, period = this.sequenceNumberPeriod;
        for (let i = 0, il = sw.length; i < il; ++i) {
            const label = sw.residueLabel(i);
            // add sequence number before name so the html element do not get separated by a line-break
            if (hasNumbers && i % period === 0 && i < il) {
                elems[elems.length] = this.getSequenceNumberSpan(i, label);
            }
            elems[elems.length] = this.residue(i, label, sw.markerArray[i]);
            if (!this.props.hideSecondaryStructure) {
                const span = this.getSequenceSecondaryStructureSpan(i);
                if (!span) continue;
                elems[elems.length] = span;
            }
        }

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
