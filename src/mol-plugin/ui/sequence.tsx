/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react'
import { Structure, StructureSequence, Queries, StructureSelection, StructureProperties as SP, StructureQuery, StructureElement, Unit } from '../../mol-model/structure';
import { PluginUIComponent } from './base';
import { StateTreeSpine } from '../../mol-state/tree/spine';
import { PluginStateObject as SO } from '../state/objects';
import { Interaction } from '../util/interaction';
import { OrderedSet, Interval } from '../../mol-data/int';
import { Loci } from '../../mol-model/loci';
import { applyMarkerAction, MarkerAction } from '../../mol-geo/geometry/marker-data';
import { ButtonsType, ModifiersKeys, getButtons, getModifiers } from '../../mol-util/input/input-observer';


export class SequenceView extends PluginUIComponent<{ }, { }> {
    private spine: StateTreeSpine.Impl

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

    private getStructure() {
        const so = this.spine && this.spine.getRootOfType(SO.Molecule.Structure)
        return so && so.data
    }

    render() {
        const s = this.getStructure();
        if (!s) return <div className='msp-sequence'>
            <div className='msp-sequence-entity'>No structure available</div>
        </div>;

        const seqs = s.models[0].sequence.sequences;
        return <div className='msp-sequence'>
            {seqs.map((seq, i) => <EntitySequence key={i} seq={seq} structure={s} /> )}
        </div>;
    }
}

function createQuery(entityId: string, label_seq_id: number) {
    return Queries.generators.atoms({
        entityTest: ctx => {
            return SP.entity.id(ctx.element) === entityId
        },
        residueTest: ctx => {
            if (ctx.element.unit.kind === Unit.Kind.Atomic) {
                return SP.residue.label_seq_id(ctx.element) === label_seq_id
            } else {
                return (
                    SP.coarse.seq_id_begin(ctx.element) <= label_seq_id &&
                    SP.coarse.seq_id_end(ctx.element) >= label_seq_id
                )
            }
        }
    });
}

function getSeqIdInterval(location: StructureElement): Interval {
    const { unit, element } = location
    const { model } = unit
    switch (unit.kind) {
        case Unit.Kind.Atomic:
            const residueIndex = model.atomicHierarchy.residueAtomSegments.index[element]
            const seqId = model.atomicHierarchy.residues.label_seq_id.value(residueIndex)
            return Interval.ofSingleton(seqId)
        case Unit.Kind.Spheres:
            return Interval.ofRange(
                model.coarseHierarchy.spheres.seq_id_begin.value(element),
                model.coarseHierarchy.spheres.seq_id_end.value(element)
            )
        case Unit.Kind.Gaussians:
            return Interval.ofRange(
                model.coarseHierarchy.gaussians.seq_id_begin.value(element),
                model.coarseHierarchy.gaussians.seq_id_end.value(element)
            )
    }
}

type StructureSeq = { structure: Structure, seq: StructureSequence.Entity }

function eachResidue(loci: Loci, structureSeq: StructureSeq, apply: (interval: Interval) => boolean) {
    let changed = false
    const { structure, seq } = structureSeq
    if (!StructureElement.isLoci(loci)) return false
    if (!Structure.areParentsEquivalent(loci.structure, structure)) return false
    const l = StructureElement.create()
    for (const e of loci.elements) {
        l.unit = e.unit
        OrderedSet.forEach(e.indices, v => {
            l.element = e.unit.elements[v]
            const entityId = SP.entity.id(l)
            if (entityId === seq.entityId) {
                if (apply(getSeqIdInterval(l))) changed = true
            }
        })
    }
    return changed
}

function markResidue(loci: Loci, structureSeq: StructureSeq, array: Uint8Array, action: MarkerAction) {
    const { structure, seq } = structureSeq
    return eachResidue(loci, { structure , seq }, (i: Interval) => {
        let changed = false
        OrderedSet.forEach(i, (v: number) => {
            const start = Interval.start(i) - 1
            const end = Interval.end(i) - 1
            if (applyMarkerAction(array, start, end, action)) changed = true
        })
        return changed
    })
}

// TODO: this is really inefficient and should be done using a canvas.
class EntitySequence extends PluginUIComponent<{ seq: StructureSequence.Entity, structure: Structure }, { markerData: { array: Uint8Array } }> {
    state = {
        markerData: { array: new Uint8Array(this.props.seq.sequence.sequence.length) }
    }

    private lociHighlightProvider = (loci: Interaction.Loci, action: MarkerAction) => {
        const { array } = this.state.markerData;
        const { structure, seq } = this.props
        const changed = markResidue(loci.loci, { structure , seq }, array, action)
        if (changed) this.setState({ markerData: { array } })
    }

    private lociSelectionProvider = (loci: Interaction.Loci, action: MarkerAction) => {
        const { array } = this.state.markerData;
        const { structure, seq } = this.props
        const changed = markResidue(loci.loci, { structure , seq }, array, action)
        if (changed) this.setState({ markerData: { array } })
    }

    componentDidMount() {
        this.plugin.lociHighlights.addProvider(this.lociHighlightProvider)
        this.plugin.lociSelections.addProvider(this.lociSelectionProvider)
    }

    componentWillUnmount() {
        this.plugin.lociHighlights.removeProvider(this.lociHighlightProvider)
        this.plugin.lociSelections.removeProvider(this.lociSelectionProvider)
    }

    getLoci(seqId: number) {
        const query = createQuery(this.props.seq.entityId, seqId);
        return StructureSelection.toLoci2(StructureQuery.run(query, this.props.structure));
    }

    highlight(seqId?: number, modifiers?: ModifiersKeys) {
        const ev = { current: Interaction.Loci.Empty, modifiers }
        if (seqId !== undefined) {
            const loci = this.getLoci(seqId);
            if (loci.elements.length > 0) ev.current = { loci };
        }
        this.plugin.behaviors.interaction.highlight.next(ev)
    }

    click(seqId: number | undefined, buttons: ButtonsType, modifiers: ModifiersKeys) {
        const ev = { current: Interaction.Loci.Empty, buttons, modifiers }
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
        const { seq } = this.props;
        const { offset, sequence } = seq.sequence;

        const elems: JSX.Element[] = [];
        for (let i = 0, _i = sequence.length; i < _i; i++) {
            elems[elems.length] = <Residue seqId={offset + i + 1} letter={sequence[i]} parent={this} marker={markerData.array[i]} key={i} />;
        }

        return <div
            className='msp-sequence-entity'
            onContextMenu={this.contextMenu}
            onMouseDown={this.mouseDown}
        >
            <span style={{ fontWeight: 'bold' }}>{this.props.seq.entityId}:{offset}&nbsp;</span>
            {elems}
        </div>;
    }
}

class Residue extends PluginUIComponent<{ seqId: number, letter: string, parent: EntitySequence, marker: number }> {

    mouseEnter = (e: React.MouseEvent) => {
        const modifiers = getModifiers(e.nativeEvent)
        this.props.parent.highlight(this.props.seqId, modifiers);
    }

    mouseLeave = () => {
        this.props.parent.highlight();
    }

    mouseDown = (e: React.MouseEvent) => {
        const buttons = getButtons(e.nativeEvent)
        const modifiers = getModifiers(e.nativeEvent)
        this.props.parent.click(this.props.seqId, buttons, modifiers);
        e.stopPropagation() // so that `parent.mouseDown` is not called
    }

    getBackgroundColor() {
        // TODO make marker color configurable
        if (this.props.marker === 0) return ''
        if (this.props.marker % 2 === 0) return 'rgb(51, 255, 25)' // selected
        return 'rgb(255, 102, 153)' // highlighted
    }

    render() {
        return <span
            onMouseEnter={this.mouseEnter}
            onMouseLeave={this.mouseLeave}
            onMouseDown={this.mouseDown}
            style={{ backgroundColor: this.getBackgroundColor() }}>
            {this.props.letter}
        </span>;
    }
}