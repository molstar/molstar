/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { shallowMergeArray } from '../mol-util/object';
import { RxEventHelper } from '../mol-util/rx-event-helper';

export class PluginComponent<State> {
    private _ev: RxEventHelper | undefined;

    protected get ev() {
        return this._ev || (this._ev = RxEventHelper.create());
    }

    private _state: State;

    protected updateState(...states: Partial<State>[]): boolean {
        const latest = this.state;
        const s = shallowMergeArray(latest, states);
        if (s !== latest) {
            this._state = s;
            return true;
        }
        return false;
    }

    get state() {
        return this._state;
    }

    dispose() {
        if (this._ev) this._ev.dispose();
    }

    constructor(initialState: State) {
        this._state = initialState;
    }
}