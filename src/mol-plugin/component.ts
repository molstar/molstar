/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BehaviorSubject, Observable, Subject } from 'rxjs';
import { PluginContext } from './context';
import { shallowMergeArray } from 'mol-util/object';

export class PluginComponent<State> {
    private _state: BehaviorSubject<State>;
    private _updated = new Subject();

    updateState(...states: Partial<State>[]): boolean {
        const latest = this.latestState;
        const s = shallowMergeArray(latest, states);
        if (s !== latest) {
            this._state.next(s);
            return true;
        }
        return false;
    }

    get states() {
        return <Observable<State>>this._state;
    }

    get latestState() {
        return this._state.value;
    }

    get updated() {
        return <Observable<{}>>this._updated;
    }

    triggerUpdate() {
        this._updated.next({});
    }

    constructor(public context: PluginContext, initialState: State) {
        this._state = new BehaviorSubject<State>(initialState);
    }
}