/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react'
import { PluginUIComponent } from './base';
import { PluginStateObject as PSO } from '../state/objects';
import { Sequence } from './sequence/sequence';
import { Structure, StructureElement, StructureProperties as SP, Unit } from '../../mol-model/structure';
import { SequenceWrapper } from './sequence/wrapper';
import { PolymerSequenceWrapper } from './sequence/polymer';
import { StructureElementSelectionManager } from '../util/structure-element-selection';
import { MarkerAction } from '../../mol-util/marker-action';
import { ParameterControls } from './controls/parameters';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { HeteroSequenceWrapper } from './sequence/hetero';
import { State, StateSelection } from '../../mol-state';

function opKey(l: StructureElement.Location) {
    const ids = SP.unit.pdbx_struct_oper_list_ids(l)
    const ncs = SP.unit.struct_ncs_oper_id(l)
    const hkl = SP.unit.hkl(l)
    const spgrOp = SP.unit.spgrOp(l)
    return `${ids.sort().join(',')}|${ncs}|${hkl}|${spgrOp}`
}

function getSequenceWrapper(state: SequenceViewState, structureSelection: StructureElementSelectionManager): SequenceWrapper.Any | undefined {
    const { structure, entityId, invariantUnitId, operatorKey } = state
    const l = StructureElement.Location.create()
    for (const unit of structure.units) {
        StructureElement.Location.set(l, unit, unit.elements[0])
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
    const l = StructureElement.Location.create()
    const seen = new Set<string>()

    for (const unit of structure.units) {
        StructureElement.Location.set(l, unit, unit.elements[0])
        const id = SP.entity.id(l)
        if (seen.has(id)) continue

        const label = `${id}: ${SP.entity.pdbx_description(l).join(', ')}`
        options.push([ id, label ])
        seen.add(id)
    }

    if (options.length === 0) options.push(['', 'No entities'])
    return options
}

function getUnitOptions(structure: Structure, entityId: string) {
    const options: [number, string][] = []
    const l = StructureElement.Location.create()
    const seen = new Set<number>()
    const water = new Map<string, number>()

    for (const unit of structure.units) {
        StructureElement.Location.set(l, unit, unit.elements[0])
        if (SP.entity.id(l) !== entityId) continue

        const id = unit.invariantId
        if (seen.has(id)) continue

        // TODO handle special cases
        // - more than one chain in a unit
        // - chain spread over multiple units
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

    if (options.length === 0) options.push([-1, 'No units'])
    return options
}

function getOperatorOptions(structure: Structure, entityId: string, invariantUnitId: number) {
    const options: [string, string][] = []
    const l = StructureElement.Location.create()
    const seen = new Set<string>()

    for (const unit of structure.units) {
        StructureElement.Location.set(l, unit, unit.elements[0])
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

function getStructureOptions(state: State) {
    const options: [string, string][] = []

    const structures = state.select(StateSelection.Generators.rootsOfType(PSO.Molecule.Structure))
    for (const s of structures) {
        options.push([s.transform.ref, s.obj!.data.label])
    }

    if (options.length === 0) options.push(['', 'No structure'])
    return options
}

type SequenceViewState = {
    structure: Structure,
    structureRef: string,
    entityId: string,
    invariantUnitId: number,
    operatorKey: string
}

export class SequenceView extends PluginUIComponent<{ }, SequenceViewState> {
    state = { structure: Structure.Empty, structureRef: '', entityId: '', invariantUnitId: -1, operatorKey: '' }

    componentDidMount() {
        if (this.plugin.state.dataState.select(StateSelection.Generators.rootsOfType(PSO.Molecule.Structure)).length > 0) this.setState(this.getInitialState())

        this.subscribe(this.plugin.events.state.object.updated, ({ ref, obj }) => {
            if (ref === this.state.structureRef && obj && obj.type === PSO.Molecule.Structure.type && obj.data !== this.state.structure) {
                this.setState(this.getInitialState())
            }
        });

        this.subscribe(this.plugin.events.state.object.created, ({ obj }) => {
            if (obj && obj.type === PSO.Molecule.Structure.type) {
                this.setState(this.getInitialState())
            }
        });

        this.subscribe(this.plugin.events.state.object.removed, ({ obj }) => {
            if (obj && obj.type === PSO.Molecule.Structure.type && obj.data === this.state.structure) {
                this.setState(this.getInitialState())
            }
        });
    }

    private getStructure(ref: string) {
        const state = this.plugin.state.dataState;
        const cell = state.select(ref)[0];
        if (!ref || !cell || !cell.obj) return Structure.Empty;
        return (cell.obj as PSO.Molecule.Structure).data;
    }

    private getSequenceWrapper() {
        return getSequenceWrapper(this.state, this.plugin.helpers.structureSelectionManager)
    }

    private getInitialState(): SequenceViewState {
        // TODO reuse selected values from previous state if applicable
        const structureRef = getStructureOptions(this.plugin.state.dataState)[0][0]
        const structure = this.getStructure(structureRef)
        const entityId = getEntityOptions(structure)[0][0]
        const invariantUnitId = getUnitOptions(structure, entityId)[0][0]
        const operatorKey = getOperatorOptions(structure, entityId, invariantUnitId)[0][0]
        return { structure, structureRef, entityId, invariantUnitId, operatorKey }
    }

    private get params() {
        const { structure, entityId, invariantUnitId } = this.state
        const structureOptions = getStructureOptions(this.plugin.state.dataState)
        const entityOptions = getEntityOptions(structure)
        const unitOptions = getUnitOptions(structure, entityId)
        const operatorOptions = getOperatorOptions(structure, entityId, invariantUnitId)
        return {
            structure: PD.Select(structureOptions[0][0], structureOptions, { shortLabel: true }),
            entity: PD.Select(entityOptions[0][0], entityOptions, { shortLabel: true }),
            unit: PD.Select(unitOptions[0][0], unitOptions, { shortLabel: true, twoColumns: true }),
            operator: PD.Select(operatorOptions[0][0], operatorOptions, { shortLabel: true, twoColumns: true })
        }
    }

    private get values(): PD.Values<SequenceView['params']> {
        return {
            structure: this.state.structureRef,
            entity: this.state.entityId,
            unit: this.state.invariantUnitId,
            operator: this.state.operatorKey
        }
    }

    // TODO try to use selected option from previous state
    private setParamProps = (p: { param: PD.Base<any>, name: string, value: any }) => {
        const state = { ...this.state }
        // console.log(p.name, p.value)
        switch (p.name) {
            case 'structure':
                state.structureRef = p.value
                state.structure = this.getStructure(p.value)
                state.entityId = getEntityOptions(state.structure)[0][0]
                state.invariantUnitId = getUnitOptions(state.structure, state.entityId)[0][0]
                state.operatorKey = getOperatorOptions(state.structure, state.entityId, state.invariantUnitId)[0][0]
                break
            case 'entity':
                state.entityId = p.value
                state.invariantUnitId = getUnitOptions(state.structure, state.entityId)[0][0]
                state.operatorKey = getOperatorOptions(state.structure, state.entityId, state.invariantUnitId)[0][0]
                break
            case 'unit':
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
        if (this.getStructure(this.state.structureRef) === Structure.Empty) return <div className='msp-sequence'>
            <div className='msp-sequence-wrapper'>No structure available</div>
        </div>;

        const sequenceWrapper = this.getSequenceWrapper()

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