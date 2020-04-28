/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { shallowMergeArray } from '../mol-util/object';
import { RxEventHelper } from '../mol-util/rx-event-helper';
import { Subscription, Observable } from 'rxjs';
import { arraySetRemove } from '../mol-util/array';

export class PluginComponent {
    private _ev: RxEventHelper | undefined;
    private subs: Subscription[] | undefined = void 0;

    protected subscribe<T>(obs: Observable<T>, action: (v: T) => void) {
        if (typeof this.subs === 'undefined') this.subs = [];

        let sub: Subscription | undefined = obs.subscribe(action);
        this.subs.push(sub);

        return {
            unsubscribe: () => {
                if (sub && this.subs && arraySetRemove(this.subs, sub)) {
                    sub.unsubscribe();
                    sub = void 0;
                }
            }
        };
    }

    protected get ev() {
        return this._ev || (this._ev = RxEventHelper.create());
    }

    dispose() {
        if (this._ev) this._ev.dispose();
        if (this.subs) {
            for (const s of this.subs) s.unsubscribe();
            this.subs = void 0;
        }
    }
}

export class StatefulPluginComponent<State> extends PluginComponent {
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

    constructor(initialState: State) {
        super();
        this._state = initialState;
    }
}