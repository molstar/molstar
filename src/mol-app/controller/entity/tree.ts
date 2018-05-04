/*
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol
 * Copyright (c) 2016 - now David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */

import { Context } from '../../context/context'
import { Controller } from '../controller';
import { AnyEntity } from 'mol-view/state/entity';

export interface EntityTreeState {
    entities: Set<AnyEntity>
}

export class EntityTreeController extends Controller<EntityTreeState> {
    constructor(context: Context) {
        super(context, { entities: new Set() });

        context.stage.ctx.change.subscribe(() => {
            if (context.stage.ctx) {
                this.state.next({ entities: context.stage.ctx.entities }) // TODO
                this.setState({ entities: context.stage.ctx.entities })
            }
        })
    }
}