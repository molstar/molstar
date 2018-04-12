/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { Observable, Subscription } from 'rxjs';

export default class Observer<S, P> extends React.Component<S, P> {
    private _subs: Subscription[] = []

    subscribe<T>(obs: Observable<T>, onNext: (v: T) => void) {
        this._subs.push(obs.subscribe(onNext));
    }

    componentWillUnmount() {
        for (const s of this._subs) s.unsubscribe();
        this._subs = [];
    }
}