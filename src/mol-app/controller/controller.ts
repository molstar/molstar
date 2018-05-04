/*
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol
 * Copyright (c) 2016 - now David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */

import { BehaviorSubject } from 'rxjs';
import { merge } from 'mol-util';
import { Context } from '../context/context'
import { LayoutRegion } from './layout';

export class Controller<State> {

    private _state = new BehaviorSubject<State>(<any>void 0);
    private _latestState: State = <any>void 0;

    get dispatcher() {
        return this.context.dispatcher;
    }

    setState(...states: Partial<State>[]) {
        let s = merge(this._latestState, ...states);
        if (s !== this._latestState) {
            this._latestState = s;
            this._state.next(s);
        }
    }

    get state() {
        return this._state;
    }

    get latestState() {
        return this._latestState;
    }

    constructor(public context: Context, initialState: State) {
        this._latestState = initialState;
    }
}

export interface ControllerInfo {
    key: string;
    controller: Controller<any>;
    view: any;
    region: LayoutRegion;
    isStatic?: boolean;
}