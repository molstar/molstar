/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { Observable, Subscription } from 'rxjs';
import { PluginContext } from '../context';

export const PluginReactContext = React.createContext(void 0 as any as PluginContext);

export abstract class PluginComponent<P = {}, S = {}, SS = {}> extends React.Component<P, S, SS> {
    static contextType = PluginReactContext;
    readonly plugin: PluginContext;

    private subs: Subscription[] | undefined = void 0;

    protected subscribe<T>(obs: Observable<T>, action: (v: T) => void) {
        if (typeof this.subs === 'undefined') this.subs = []
        this.subs.push(obs.subscribe(action));
    }

    componentWillUnmount() {
        if (!this.subs) return;
        for (const s of this.subs) s.unsubscribe();
    }

    protected init?(): void;

    constructor(props: P, context?: any) {
        super(props, context);
        this.plugin = context;
        if (this.init) this.init();
    }
}

export abstract class PurePluginComponent<P = {}, S = {}, SS = {}> extends React.PureComponent<P, S, SS> {
    static contextType = PluginReactContext;
    readonly plugin: PluginContext;

    private subs: Subscription[] | undefined = void 0;

    protected subscribe<T>(obs: Observable<T>, action: (v: T) => void) {
        if (typeof this.subs === 'undefined') this.subs = []
        this.subs.push(obs.subscribe(action));
    }

    componentWillUnmount() {
        if (!this.subs) return;
        for (const s of this.subs) s.unsubscribe();
    }

    protected init?(): void;

    constructor(props: P, context?: any) {
        super(props, context);
        this.plugin = context;
        if (this.init) this.init();
    }
}