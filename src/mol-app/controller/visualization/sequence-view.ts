/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { shallowClone } from 'mol-util';
import { Context } from '../../context/context'
import { Controller } from '../controller';
import { Structure } from 'mol-model/structure';

export const DefaultSequenceViewState = {
    structure: void 0 as (Structure | undefined)
}
export type SequenceViewState = typeof DefaultSequenceViewState

export class SequenceViewController extends Controller<SequenceViewState> {
    constructor(context: Context) {
        super(context, shallowClone(DefaultSequenceViewState));
    }
}