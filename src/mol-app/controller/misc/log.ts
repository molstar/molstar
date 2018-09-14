/*
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol
 * Copyright (c) 2016 - now David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */

import produce from 'immer'

import { Controller } from '../controller'
import { Context } from '../../context/context';
import { LogEvent } from '../../event/basic';
import { Logger } from '../../service/logger';

export class LogController extends Controller<{ entries: Logger.Entry[] }> {
    constructor(context: Context) {
        super(context, { entries: [] });

        LogEvent.getStream(this.context)
            .subscribe(e => this.setState({
                entries: produce(this.latestState.entries, _entries => { _entries.push(e.data) })
            }))
    }
}