/*
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol
 * Copyright (c) 2016 - now David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */

import { Context } from '../../context/context'
import { Controller } from '../controller';
import { AnyTransform } from 'mol-view/state/transform';
import { AnyEntity } from 'mol-view/state/entity';

export interface TransformListState {
    entity?: AnyEntity
    transforms: AnyTransform[]
}

export class TransformListController extends Controller<TransformListState> {
    constructor(context: Context) {
        super(context, { transforms: [], entity: undefined });

        context.currentTransforms.subscribe((transforms) => {
            this.state.next({ transforms, entity: context.currentEntity.getValue() }) // TODO
            this.setState({ transforms, entity: context.currentEntity.getValue() })
        })
    }
}