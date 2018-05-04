/*
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol
 * Copyright (c) 2016 - now David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */

import * as React from 'react'
import { Observable, Subscription } from 'rxjs';
import { merge, shallowEqual } from 'mol-util'
import { Context } from '../context/context';
import { Controller } from '../controller/controller';

export abstract class PureView<State, Props, ViewState> extends React.Component<{
    state: State
    onChange: (s: State) => void
} & Props, ViewState> {

    protected update(s: State) {
        let ns = merge<State>(this.props.state, s);
        if (ns !== this.props.state as any) this.props.onChange(ns);
    }

    shouldComponentUpdate(nextProps: any, nextState: any) {
        return !shallowEqual(this.props, nextProps) || !shallowEqual(this.state, nextState);
    }
}

export abstract class ComponentView<Props> extends React.Component<{ context: Context } & Props, {}> {

    // shouldComponentUpdate(nextProps: any, nextState: any) {
    //     return !shallowEqual(this.props, nextProps);
    // }

    private subs: Subscription[] = [];
    protected subscribe<T>(stream: Observable<T>, obs: (n: T) => void) {
        let sub = stream.subscribe(obs);
        this.subs.push(sub);
        return sub;
    }

    protected unsubscribe(sub: Subscription) {
        let idx = this.subs.indexOf(sub);
        for (let i = idx; i < this.subs.length - 1; i++) {
            this.subs[i] = this.subs[i + 1];
        }
        sub.unsubscribe();
        this.subs.pop();
    }

    componentWillUnmount() {
        for (let s of this.subs) s.unsubscribe();
        this.subs = [];
    }
}

export abstract class ObserverView<P, S> extends React.Component<P, S> {
    private subs: Subscription[] = [];

    protected subscribe<T>(stream: Observable<T>, obs: (n: T) => void) {
        let sub = stream.subscribe(obs);
        this.subs.push(sub);
        return sub;
    }

    protected unsubscribe(sub: Subscription) {
        let idx = this.subs.indexOf(sub);
        for (let i = idx; i < this.subs.length - 1; i++) {
            this.subs[i] = this.subs[i + 1];
        }
        sub.unsubscribe();
        this.subs.pop();
    }

    componentWillUnmount() {
        for (let s of this.subs) s.unsubscribe();
        this.subs = [];
    }
}

export abstract class View<T extends Controller<any>, State, CustomProps>
    extends ObserverView<{ controller: T } & CustomProps, State> {

    public get controller(): T {
        return this.props.controller as any;
    }

    componentWillMount() {
        this.subscribe(this.controller.state as any, (s) => {
            this.forceUpdate()
        });
    }
}