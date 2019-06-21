/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react'
import { PluginUIComponent } from './base';
import { StateTreeSpine } from '../../mol-state/tree/spine';
import { PluginStateObject as SO } from '../state/objects';
import { Sequence } from './sequence/sequence';
import { Structure, StructureElement, StructureProperties as SP } from '../../mol-model/structure';
import { SequenceWrapper } from './sequence/util';
import { PolymerSequenceWrapper } from './sequence/polymer';
import { StructureElementSelectionManager } from '../util/structure-element-selection';
import { MarkerAction } from '../../mol-util/marker-action';
import { ParameterControls } from './controls/parameters';
import { ParamDefinition as PD } from '../../mol-util/param-definition';

function opKey(l: StructureElement) {
    const ids = SP.unit.pdbx_struct_oper_list_ids(l)
    const ncs = SP.unit.struct_ncs_oper_id(l)
    const hkl = SP.unit.hkl(l)
    const spgrOp = SP.unit.spgrOp(l)
    return `${ids.sort().join(',')}|${ncs}|${hkl}|${spgrOp}`
}

function getSequenceWrapper(state: SequenceViewState, structureSelection: StructureElementSelectionManager): SequenceWrapper.Any | undefined {
    const { structure, entity, chain, operator } = state
    const l = StructureElement.create()
    for (let i = 0, il = structure.units.length; i < il; ++i) {
        const unit = structure.units[i]
        if (unit.polymerElements.length === 0) continue

        StructureElement.set(l, unit, unit.elements[0])
        if (SP.entity.id(l) !== entity) continue
        if (SP.chain.label_asym_id(l) !== chain) continue
        if (opKey(l) !== operator) continue

        // console.log('new PolymerSequenceWrapper', structureSelection.get(structure))
        const sw = new PolymerSequenceWrapper({ structure, unit })
        sw.markResidue(structureSelection.get(structure), MarkerAction.Select)
        return sw
    }
}

function getEntityOptions(structure: Structure) {
    const options: [string, string][] = []
    const l = StructureElement.create()
    const seen = new Set<string>()

    structure.units.forEach(unit => {
        if (unit.polymerElements.length === 0) return

        StructureElement.set(l, unit, unit.elements[0])
        const id = SP.entity.id(l)
        if (seen.has(id)) return

        const label = `${id}: ${SP.entity.pdbx_description(l).join(', ')}`
        options.push([ id, label ])
        seen.add(id)
    })

    if (options.length === 0) options.push(['', 'No entities'])
    return options
}

function getChainOptions(structure: Structure, entityId: string) {
    const options: [string, string][] = []
    const l = StructureElement.create()
    const seen = new Set<string>()

    structure.units.forEach(unit => {
        if (unit.polymerElements.length === 0) return

        StructureElement.set(l, unit, unit.elements[0])
        if (SP.entity.id(l) !== entityId) return

        const id = SP.chain.label_asym_id(l)
        if (seen.has(id)) return

        const label = `${id}: ${SP.chain.auth_asym_id(l)}`
        options.push([ id, label ])
        seen.add(id)
    })

    if (options.length === 0) options.push(['', 'No chains'])
    return options
}

function getOperatorOptions(structure: Structure, entityId: string, label_asym_id: string) {
    const options: [string, string][] = []
    const l = StructureElement.create()
    const seen = new Set<string>()

    structure.units.forEach(unit => {
        if (unit.polymerElements.length === 0) return
        StructureElement.set(l, unit, unit.elements[0])
        if (SP.entity.id(l) !== entityId) return
        if (SP.chain.label_asym_id(l) !== label_asym_id) return

        const id = opKey(l)
        if (seen.has(id)) return

        const label = unit.conformation.operator.name
        options.push([ id, label ])
        seen.add(id)
    })

    if (options.length === 0) options.push(['', 'No operators'])
    return options
}

type SequenceViewState = { structure: Structure, entity: string, chain: string, operator: string }

export class SequenceView extends PluginUIComponent<{ }, SequenceViewState> {
    private spine: StateTreeSpine.Impl

    state = { structure: Structure.Empty, entity: '', chain: '', operator: '' }

    constructor(props: {}, context?: any) {
        super(props, context);
        this.spine = new StateTreeSpine.Impl(this.plugin.state.dataState.cells);
    }

    componentDidMount() {
        this.subscribe(this.plugin.state.behavior.currentObject, o => {
            const current = this.plugin.state.dataState.cells.get(o.ref)!;
            this.spine.current = current
            this.setState(this.getInitialState())
        });

        this.subscribe(this.plugin.events.state.object.updated, ({ ref, state }) => {
            const current = this.spine.current;
            if (!current || current.sourceRef !== ref || current.state !== state) return;
            this.setState(this.getInitialState())
        });
    }

    private getStructure() {
        const so = this.spine.getRootOfType(SO.Molecule.Structure)
        return (so && so.data) || Structure.Empty
    }

    private getSequenceWrapper() {
        return getSequenceWrapper(this.state, this.plugin.helpers.structureSelection)
    }

    private getInitialState(): SequenceViewState {
        const structure = this.getStructure()
        const entity = getEntityOptions(structure)[0][0]
        const chain = getChainOptions(structure, entity)[0][0]
        const operator = getOperatorOptions(structure, entity, chain)[0][0]
        return { structure, entity, chain, operator }
    }

    private get params() {
        const { structure, entity, chain } = this.state
        const entityOptions = getEntityOptions(structure)
        const chainOptions = getChainOptions(structure, entity)
        const operatorOptions = getOperatorOptions(structure, entity, chain)
        return {
            entity: PD.Select(entityOptions[0][0], entityOptions),
            chain: PD.Select(chainOptions[0][0], chainOptions),
            operator: PD.Select(operatorOptions[0][0], operatorOptions)
        }
    }

    private setParamProps = (p: { param: PD.Base<any>, name: string, value: any }) => {
        const state = { ...this.state }
        switch (p.name) {
            case 'entity':
                state.entity = p.value
                state.chain = getChainOptions(state.structure, state.entity)[0][0]
                state.operator = getOperatorOptions(state.structure, state.entity, state.chain)[0][0]
                break
            case 'chain':
                state.chain = p.value
                state.operator = getOperatorOptions(state.structure, state.entity, state.chain)[0][0]
                break
            case 'operator':
                state.operator = p.value
                break
        }
        this.setState(state)
    }

    render() {
        if (this.state.structure === Structure.Empty) return <div className='msp-sequence'>
            <div className='msp-sequence-wrapper'>No structure available</div>
        </div>;

        const sequenceWrapper = this.getSequenceWrapper()
        return <div className='msp-sequence'>
            <div className='msp-sequence-select'>
                <ParameterControls params={this.params} values={this.state} onChange={this.setParamProps} />
            </div>
            {sequenceWrapper !== undefined
                ? <Sequence sequenceWrapper={sequenceWrapper} />
                : <div className='msp-sequence-wrapper'>No sequence available</div>}
        </div>;
    }
}