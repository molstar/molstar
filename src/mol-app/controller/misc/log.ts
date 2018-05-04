/*
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol
 * Copyright (c) 2016 - now David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */

import { List } from 'immutable'

import { Controller } from '../controller'
import { Context } from '../../context/context';
import { LogEvent } from '../../event/basic';
import { Logger } from '../../service/logger';

export class LogController extends Controller<{ entries: List<Logger.Entry> }> {
    constructor(context: Context) {
        super(context, { entries: List<Logger.Entry>() });

        LogEvent.getStream(this.context)
            .subscribe(e => this.setState({ entries: this.latestState.entries.push(e.data) }))
    }
}