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
import { Structure, StructureElement, StructureProperties as SP, Unit } from '../../mol-model/structure';
import { SequenceWrapper } from './sequence/wrapper';
import { PolymerSequenceWrapper } from './sequence/polymer';
import { StructureElementSelectionManager } from '../util/structure-element-selection';
import { MarkerAction } from '../../mol-util/marker-action';
import { ParameterControls } from './controls/parameters';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { HeteroSequenceWrapper } from './sequence/hetero';

function opKey(l: StructureElement) {
    const ids = SP.unit.pdbx_struct_oper_list_ids(l)
    const ncs = SP.unit.struct_ncs_oper_id(l)
    const hkl = SP.unit.hkl(l)
    const spgrOp = SP.unit.spgrOp(l)
    return `${ids.sort().join(',')}|${ncs}|${hkl}|${spgrOp}`
}

function getSequenceWrapper(state: SequenceViewState, structureSelection: StructureElementSelectionManager): SequenceWrapper.Any | undefined {
    const { structure, entityId, invariantUnitId, operatorKey } = state
    const l = StructureElement.create()
    for (const unit of structure.units) {
        StructureElement.set(l, unit, unit.elements[0])
        if (SP.entity.id(l) !== entityId) continue
        if (unit.invariantId !== invariantUnitId) continue
        if (opKey(l) !== operatorKey) continue

        const Wrapper = unit.polymerElements.length ? PolymerSequenceWrapper : HeteroSequenceWrapper
        const sw = new Wrapper({ structure, unit })
        sw.markResidue(structureSelection.get(structure), MarkerAction.Select)
        return sw
    }
}

function getEntityOptions(structure: Structure) {
    const options: [string, string][] = []
    const l = StructureElement.create()
    const seen = new Set<string>()

    for (const unit of structure.units) {
        StructureElement.set(l, unit, unit.elements[0])
        const id = SP.entity.id(l)
        if (seen.has(id)) continue

        const label = `${id}: ${SP.entity.pdbx_description(l).join(', ')}`
        options.push([ id, label ])
        seen.add(id)
    }

    if (options.length === 0) options.push(['', 'No entities'])
    return options
}

function getChainOptions(structure: Structure, entityId: string) {
    const options: [number, string][] = []
    const l = StructureElement.create()
    const seen = new Set<number>()
    const water = new Map<string, number>()

    for (const unit of structure.units) {
        StructureElement.set(l, unit, unit.elements[0])
        if (SP.entity.id(l) !== entityId) continue

        const id = unit.invariantId
        if (seen.has(id)) continue

        let label = ''
        if (Unit.isAtomic(unit)) {
            label = `${SP.chain.label_asym_id(l)}: ${SP.chain.auth_asym_id(l)}`
        } else {
            label = `${SP.coarse.asym_id(l)}`
        }
        if (SP.entity.type(l) === 'water') {
            const count = water.get(label) || 1
            water.set(label, count + 1)
            label += ` #${count}`
        }

        options.push([ id, label ])
        seen.add(id)
    }

    if (options.length === 0) options.push([-1, 'No chains'])
    return options
}

function getOperatorOptions(structure: Structure, entityId: string, invariantUnitId: number) {
    const options: [string, string][] = []
    const l = StructureElement.create()
    const seen = new Set<string>()

    for (const unit of structure.units) {
        StructureElement.set(l, unit, unit.elements[0])
        if (SP.entity.id(l) !== entityId) continue
        if (unit.invariantId !== invariantUnitId) continue

        const id = opKey(l)
        if (seen.has(id)) continue

        const label = unit.conformation.operator.name
        options.push([ id, label ])
        seen.add(id)
    }

    if (options.length === 0) options.push(['', 'No operators'])
    return options
}

type SequenceViewState = { structure: Structure, entityId: string, invariantUnitId: number, operatorKey: string }

export class SequenceView extends PluginUIComponent<{ }, SequenceViewState> {
    private spine: StateTreeSpine.Impl

    state = { structure: Structure.Empty, entityId: '', invariantUnitId: -1, operatorKey: '' }

    constructor(props: {}, context?: any) {
        super(props, context);
        this.spine = new StateTreeSpine.Impl(this.plugin.state.dataState.cells);
    }

    componentDidMount() {
        this.subscribe(this.plugin.state.behavior.currentObject, o => {
            const current = this.plugin.state.dataState.cells.get(o.ref)!;
            this.spine.current = current
            if (!Structure.areParentsEqual(this.state.structure, this.getStructure())) {
                this.setState(this.getInitialState())
            }
        });

        this.subscribe(this.plugin.events.state.object.updated, ({ ref, state }) => {
            const current = this.spine.current;
            if (!current || current.sourceRef !== ref) return;
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
        const entityId = getEntityOptions(structure)[0][0]
        const invariantUnitId = getChainOptions(structure, entityId)[0][0]
        const operatorKey = getOperatorOptions(structure, entityId, invariantUnitId)[0][0]
        return { structure, entityId, invariantUnitId, operatorKey }
    }

    private get params() {
        const { structure, entityId, invariantUnitId } = this.state
        const entityOptions = getEntityOptions(structure)
        const chainOptions = getChainOptions(structure, entityId)
        const operatorOptions = getOperatorOptions(structure, entityId, invariantUnitId)
        return {
            entity: PD.Select(entityOptions[0][0], entityOptions),
            chain: PD.Select(chainOptions[0][0], chainOptions),
            operator: PD.Select(operatorOptions[0][0], operatorOptions)
        }
    }

    private get values(): PD.Values<SequenceView['params']> {
        return {
            entity: this.state.entityId,
            chain: this.state.invariantUnitId,
            operator: this.state.operatorKey
        }
    }

    // TODO try to use selected option from previous state
    private setParamProps = (p: { param: PD.Base<any>, name: string, value: any }) => {
        const state = { ...this.state }
        switch (p.name) {
            case 'entity':
                state.entityId = p.value
                state.invariantUnitId = getChainOptions(state.structure, state.entityId)[0][0]
                state.operatorKey = getOperatorOptions(state.structure, state.entityId, state.invariantUnitId)[0][0]
                break
            case 'chain':
                state.invariantUnitId = p.value
                state.operatorKey = getOperatorOptions(state.structure, state.entityId, state.invariantUnitId)[0][0]
                break
            case 'operator':
                state.operatorKey = p.value
                break
        }
        this.setState(state)
    }

    render() {
        if (this.state.structure === Structure.Empty) return <div className='msp-sequence'>
            <div className='msp-sequence-wrapper'>No structure available</div>
        </div>;

        const sequenceWrapper = this.getSequenceWrapper()

        sequenceWrapper

        return <div className='msp-sequence'>
            <div className='msp-sequence-select'>
                <ParameterControls params={this.params} values={this.values} onChange={this.setParamProps} />
            </div>
            {sequenceWrapper !== undefined
                ? (sequenceWrapper.length <= 10000
                    ? <Sequence sequenceWrapper={sequenceWrapper} />
                    : <div className='msp-sequence-wrapper'>Sequence too long</div>
                )
                : <div className='msp-sequence-wrapper'>No sequence available</div>}
        </div>;
    }
}