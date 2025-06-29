/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import * as React from 'react';
import { Subject } from 'rxjs';
import { throttleTime } from 'rxjs/operators';
import { OrderedSet } from '../../mol-data/int';
import { EveryLoci } from '../../mol-model/loci';
import { StructureElement, StructureProperties, Unit } from '../../mol-model/structure';
import { PluginCommands } from '../../mol-plugin/commands';
import { Representation } from '../../mol-repr/representation';
import { Color } from '../../mol-util/color';
import { ButtonsType, getButton, getButtons, getModifiers, ModifiersKeys } from '../../mol-util/input/input-observer';
import { MarkerAction } from '../../mol-util/marker-action';
import { PluginUIComponent } from '../base';
import { SequenceWrapper } from './wrapper';


type SequenceProps = {
    sequenceWrapper: SequenceWrapper.Any,
    sequenceNumberPeriod?: number,
    hideSequenceNumbers?: boolean,
}

/** Note, if this is changed, the CSS for `msp-sequence-number` needs adjustment too */
const MaxSequenceNumberSize = 5;

const DefaultMarkerColors = {
    selected: 'rgb(51, 255, 25)',
    highlighted: 'rgb(255, 102, 153)',
    focused: '',
};

// TODO: this is somewhat inefficient and should be done using a canvas.
export class Sequence<P extends SequenceProps> extends PluginUIComponent<P> {
    protected parentDiv = React.createRef<HTMLDivElement>();
    protected lastMouseOverSeqIdx = -1;
    protected highlightQueue = new Subject<{ seqIdx: number, buttons: number, button: number, modifiers: ModifiersKeys }>();
    protected markerColors = { ...DefaultMarkerColors };

    protected lociHighlightProvider = (loci: Representation.Loci, action: MarkerAction) => {
        const changed = this.props.sequenceWrapper.markResidue(loci.loci, action);
        if (changed) this.updateMarker();
    };

    protected lociSelectionProvider = (loci: Representation.Loci, action: MarkerAction) => {
        const changed = this.props.sequenceWrapper.markResidue(loci.loci, action);
        if (changed) this.updateMarker();
    };

    protected get sequenceNumberPeriod() {
        if (this.props.sequenceNumberPeriod !== undefined) {
            return this.props.sequenceNumberPeriod as number;
        }
        if (this.props.sequenceWrapper.length > 10) return 10;
        const lastSeqNum = this.getSequenceNumber(this.props.sequenceWrapper.length - 1);
        if (lastSeqNum.length > 1) return 5;
        return 1;
    }

    componentDidMount() {
        this.plugin.managers.interactivity.lociHighlights.addProvider(this.lociHighlightProvider);
        this.plugin.managers.interactivity.lociSelects.addProvider(this.lociSelectionProvider);

        this.subscribe(this.highlightQueue.pipe(throttleTime(3 * 16.666, void 0, { leading: true, trailing: true })), (e) => {
            const loci = this.getLoci(e.seqIdx < 0 ? void 0 : e.seqIdx);
            this.hover(loci, e.buttons, e.button, e.modifiers);
        });
        this.subscribe(this.plugin.managers.structure.focus.behaviors.current, focus => {
            this.updateFocus(focus?.loci);
            this.updateMarker();
        });

        this.updateColors();
        PluginCommands.Canvas3D.SetSettings.subscribe(this.plugin, () => {
            this.updateColors();
            this.updateMarker();
        });
    }

    updateColors() {
        if (this.plugin.canvas3d) {
            this.markerColors.highlighted = Color.toHexStyle(this.plugin.canvas3d.props.renderer.highlightColor);
            this.markerColors.selected = Color.toHexStyle(this.plugin.canvas3d.props.renderer.selectColor);
        } else {
            this.markerColors.highlighted = DefaultMarkerColors.highlighted;
            this.markerColors.selected = DefaultMarkerColors.selected;
        }
    }

    updateFocus(loci: StructureElement.Loci | undefined) {
        this.props.sequenceWrapper.markResidue(EveryLoci, 'unfocus');
        if (loci) {
            this.props.sequenceWrapper.markResidue(loci, 'focus');
        }
    }

    componentWillUnmount() {
        super.componentWillUnmount();

        this.plugin.managers.interactivity.lociHighlights.removeProvider(this.lociHighlightProvider);
        this.plugin.managers.interactivity.lociSelects.removeProvider(this.lociSelectionProvider);
    }

    getLoci(seqIdx: number | undefined) {
        if (seqIdx !== undefined) {
            const loci = this.props.sequenceWrapper.getLoci(seqIdx);
            if (!StructureElement.Loci.isEmpty(loci)) return loci;
        }
    }

    getSeqIdx(e: React.MouseEvent) {
        let seqIdx: number | undefined = undefined;
        const el = e.target as HTMLElement;
        if (el && el.getAttribute) {
            seqIdx = el.hasAttribute('data-seqid') ? +el.getAttribute('data-seqid')! : undefined;
        }
        return seqIdx;
    }

    hover(loci: StructureElement.Loci | undefined, buttons: ButtonsType, button: ButtonsType.Flag, modifiers: ModifiersKeys) {
        const ev = { current: Representation.Loci.Empty, buttons, button, modifiers };
        if (loci !== undefined && !StructureElement.Loci.isEmpty(loci)) {
            ev.current = { loci };
            if (this.mouseDownLoci) {
                const ref = this.mouseDownLoci.elements[0];
                const ext = loci.elements[0];
                const min = Math.min(OrderedSet.min(ref.indices), OrderedSet.min(ext.indices));
                const max = Math.max(OrderedSet.max(ref.indices), OrderedSet.max(ext.indices));

                const range = StructureElement.Loci(loci.structure, [{
                    unit: ref.unit,
                    indices: OrderedSet.ofRange(min as StructureElement.UnitIndex, max as StructureElement.UnitIndex)
                }]);
                ev.current = { loci: range };
            }
        }
        this.plugin.behaviors.interaction.hover.next(ev);
    }

    click(loci: StructureElement.Loci | undefined, buttons: ButtonsType, button: ButtonsType.Flag, modifiers: ModifiersKeys) {
        const ev = { current: Representation.Loci.Empty, buttons, button, modifiers };
        if (loci !== undefined && !StructureElement.Loci.isEmpty(loci)) {
            ev.current = { loci };
        }
        this.plugin.behaviors.interaction.click.next(ev);
    }

    contextMenu = (e: React.MouseEvent) => {
        e.preventDefault();
    };

    protected mouseDownLoci: StructureElement.Loci | undefined = undefined;

    mouseDown = (e: React.MouseEvent) => {
        e.stopPropagation();

        const seqIdx = this.getSeqIdx(e);
        const loci = this.getLoci(seqIdx);
        this.mouseDownLoci = loci;
    };

    mouseUp = (e: React.MouseEvent) => {
        e.stopPropagation();

        // ignore mouse-up events without a bound loci
        if (this.mouseDownLoci === undefined) return;

        const seqIdx = this.getSeqIdx(e);
        const loci = this.getLoci(seqIdx);

        if (loci) {
            const buttons = getButtons(e.nativeEvent);
            const button = getButton(e.nativeEvent);
            const modifiers = getModifiers(e.nativeEvent);

            let range = loci;
            if (!StructureElement.Loci.areEqual(this.mouseDownLoci, loci)) {
                const ref = this.mouseDownLoci.elements[0];
                const ext = loci.elements[0];
                const min = Math.min(OrderedSet.min(ref.indices), OrderedSet.min(ext.indices));
                const max = Math.max(OrderedSet.max(ref.indices), OrderedSet.max(ext.indices));

                range = StructureElement.Loci(loci.structure, [{
                    unit: ref.unit,
                    indices: OrderedSet.ofRange(min as StructureElement.UnitIndex, max as StructureElement.UnitIndex)
                }]);
            }

            this.click(range, buttons, button, modifiers);
        }
        this.mouseDownLoci = undefined;
    };

    protected getBackgroundColor(seqIdx: number) {
        const seqWrapper = this.props.sequenceWrapper;
        if (seqWrapper.isHighlighted(seqIdx)) return this.markerColors.highlighted;
        if (seqWrapper.isSelected(seqIdx)) return this.markerColors.selected;
        if (seqWrapper.isFocused(seqIdx)) return this.markerColors.focused;
        return '';
    }

    protected getResidueClass(seqIdx: number, label: string) {
        const seqWrapper = this.props.sequenceWrapper;
        const classes = [seqWrapper.residueClass(seqIdx)];
        if (label.length > 1) {
            classes.push(seqIdx === 0 ? 'msp-sequence-residue-long-begin' : 'msp-sequence-residue-long');
        }
        if (seqWrapper.isHighlighted(seqIdx)) classes.push('msp-sequence-residue-highlighted');
        if (seqWrapper.isSelected(seqIdx)) classes.push('msp-sequence-residue-selected');
        if (seqWrapper.isFocused(seqIdx)) classes.push('msp-sequence-residue-focused');
        return classes.join(' ');
    }

    protected residue(seqIdx: number, label: string) {
        return <span key={seqIdx} data-seqid={seqIdx} style={{ backgroundColor: this.getBackgroundColor(seqIdx) }} className={this.getResidueClass(seqIdx, label)}>{`\u200b${label}\u200b`}</span>;
    }

    protected getSequenceNumberClass(seqIdx: number, seqNum: string, label: string) {
        const classList = ['msp-sequence-number'];
        if (seqNum.startsWith('-')) {
            if (label.length > 1 && seqIdx > 0) classList.push('msp-sequence-number-long-negative');
            else classList.push('msp-sequence-number-negative');
        } else {
            if (label.length > 1 && seqIdx > 0) classList.push('msp-sequence-number-long');
        }
        return classList.join(' ');
    }

    protected location = StructureElement.Location.create(void 0);
    protected getSequenceNumber(seqIdx: number) {
        let seqNum = '';
        const loci = this.props.sequenceWrapper.getLoci(seqIdx);
        const l = StructureElement.Loci.getFirstLocation(loci, this.location);
        if (l) {
            if (Unit.isAtomic(l.unit)) {
                const seqId = StructureProperties.residue.auth_seq_id(l);
                const insCode = StructureProperties.residue.pdbx_PDB_ins_code(l);
                seqNum = `${seqId}${insCode ? insCode : ''}`;
            } else if (Unit.isCoarse(l.unit)) {
                seqNum = `${seqIdx + 1}`;
            }
        }
        return seqNum;
    }

    protected padSeqNum(n: string) {
        if (n.length < MaxSequenceNumberSize) return n + new Array(MaxSequenceNumberSize - n.length + 1).join('\u00A0');
        return n;
    }
    protected getSequenceNumberSpan(seqIdx: number, label: string) {
        const seqNum = this.getSequenceNumber(seqIdx);
        return <span key={`marker-${seqIdx}`} className={this.getSequenceNumberClass(seqIdx, seqNum, label)}>{this.padSeqNum(seqNum)}</span>;
    }

    protected updateMarker() {
        if (!this.parentDiv.current) return;
        const xs = this.parentDiv.current.children;
        const hasNumbers = !this.props.hideSequenceNumbers, period = this.sequenceNumberPeriod;

        const seqWrapper = this.props.sequenceWrapper;
        const seqLength = seqWrapper.length;
        let o = 0;
        for (let i = 0; i < seqLength; i++) {
            if (hasNumbers && i % period === 0 && i < seqLength) o++;
            // o + 1 to account for help icon
            const span = xs[o] as HTMLSpanElement;
            if (!span) return;
            o++;

            const className = this.getResidueClass(i, seqWrapper.residueLabel(i));
            if (span.className !== className) span.className = className;
            const backgroundColor = this.getBackgroundColor(i);
            if (span.style.backgroundColor !== backgroundColor) span.style.backgroundColor = backgroundColor;
        }
    }

    mouseMove = (e: React.MouseEvent) => {
        e.stopPropagation();

        const buttons = getButtons(e.nativeEvent);
        const button = getButton(e.nativeEvent);
        const modifiers = getModifiers(e.nativeEvent);

        const el = e.target as HTMLElement;
        if (!el || !el.getAttribute) {
            if (this.lastMouseOverSeqIdx === -1) return;
            this.lastMouseOverSeqIdx = -1;
            this.highlightQueue.next({ seqIdx: -1, buttons, button, modifiers });
            return;
        }
        const seqIdx = el.hasAttribute('data-seqid') ? +el.getAttribute('data-seqid')! : -1;
        if (this.lastMouseOverSeqIdx === seqIdx) {
            return;
        } else {
            this.lastMouseOverSeqIdx = seqIdx;
            if (this.mouseDownLoci !== undefined) {
                const loci = this.getLoci(seqIdx);
                this.hover(loci, ButtonsType.Flag.None, ButtonsType.Flag.None, modifiers);
            } else {
                this.highlightQueue.next({ seqIdx, buttons, button, modifiers });
            }
        }
    };

    mouseLeave = (e: React.MouseEvent) => {
        e.stopPropagation();
        this.mouseDownLoci = undefined;

        if (this.lastMouseOverSeqIdx === -1) return;
        this.lastMouseOverSeqIdx = -1;
        const buttons = getButtons(e.nativeEvent);
        const button = getButton(e.nativeEvent);
        const modifiers = getModifiers(e.nativeEvent);
        this.highlightQueue.next({ seqIdx: -1, buttons, button, modifiers });
    };

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
            elems[elems.length] = this.residue(i, label);
        }

        // ensure the focus markers are updated after sequenceRender is recreated
        this.updateFocus(this.plugin.managers.structure.focus.behaviors.current.value?.loci);

        // calling .updateMarker here is neccesary to ensure existing
        // residue spans are updated as react won't update them
        this.updateMarker();

        return <div
            className='msp-sequence-wrapper'
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
