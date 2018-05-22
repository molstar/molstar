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
import { Button } from '../controls/common';
import { StructureEntity } from 'mol-view/state/entity';
import { StructureCenter as StructureCenterTransform } from 'mol-view/state/transform'
import { StateContext } from 'mol-view/state/context';

export const ColorThemeInfo = {
    'atom-index': {},
    'chain-id': {},
    'element-symbol': {},
    'instance-index': {},
    'uniform': {}
}
export type ColorThemeInfo = keyof typeof ColorThemeInfo

export class StructureCenter extends View<Controller<any>, {}, { transform: StructureCenterTransform, entity: StructureEntity, ctx: StateContext }> {
    center() {
        const { ctx, entity, transform } = this.props
        transform.apply(ctx, entity)
    }

    render() {
        const { transform } = this.props

        return <div className='molstar-transformer-wrapper'>
            <div className='molstar-panel molstar-control molstar-transformer molstar-panel-expanded'>
                <div className='molstar-panel-header'>
                    <button
                        className='molstar-btn molstar-btn-link molstar-panel-expander'
                        onClick={() => {}}
                    >
                        <span>[{transform.kind}] {transform.inputKind} -> {transform.outputKind}</span>
                    </button>
                </div>
                <div className='molstar-panel-body'>
                    <div>
                        <div className='molstar-control-row molstar-options-group'>
                            <div>
                                <Button onClick={value => this.center()}>
                                    Center
                                </Button>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>;
    }
}