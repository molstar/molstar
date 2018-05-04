/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol
 * Copyright (c) 2016 - now David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'

import { View } from '../view';
import { Controller } from '../../controller/controller';
import { ModelEntity } from 'mol-view/state/entity';
import { StructureProps, ModelToStructure as ModelToStructureTransform } from 'mol-view/state/transform'
import { StateContext } from 'mol-view/state/context';

export class ModelToStructure extends View<Controller<any>, StructureProps, { transform: ModelToStructureTransform, entity: ModelEntity, ctx: StateContext }> {
    state = {
        assembly: ''
    }

    create(state?: Partial<StructureProps>) {
        const { transform, entity, ctx } = this.props
        console.log('create structure', transform, entity)
        const newState = { ...this.state, ...state }
        this.setState(newState)
        transform.apply(ctx, entity, newState)
    }

    render() {
        const { transform, entity } = this.props

        const assemblyOptions = entity.value[0].symmetry.assemblies.map((value, idx) => {
            return <option key={value.id} value={value.id}>{value.details}</option>
        })

        return <div className='molstar-transformer-wrapper'>
            <div className='molstar-panel molstar-control molstar-transformer molstar-panel-expanded'>
                <div className='molstar-panel-header'>
                    <button
                        className='molstar-btn molstar-btn-link molstar-panel-expander'
                        onClick={() => this.create()}
                    >
                        <span>[{transform.kind}] {transform.inputKind} -> {transform.outputKind}</span>
                    </button>
                </div>
                <div className='molstar-panel-body'>
                    <div>
                        <div className='molstar-control-row molstar-options-group'>
                            <span>Details</span>
                            <div>
                                <select
                                    className='molstar-form-control'
                                    value={this.state.assembly}
                                    onChange={(e) => this.create({ assembly: e.target.value })}
                                >
                                    {assemblyOptions}
                                </select>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>;
    }
}