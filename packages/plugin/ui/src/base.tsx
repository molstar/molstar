/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import type { Observable, Subscription } from 'rxjs';
import type { PluginUIContext } from '@molstar/plugin-ui/context';

export const PluginReactContext = React.createContext(void 0 as any as PluginUIContext);

export abstract class PluginUIComponent<P extends {} = {}, S = {}, SS = {}> extends React.Component<P & { children?: any }, S, SS> {
    static contextType = PluginReactContext;
    readonly plugin: PluginUIContext;

    private subs: Subscription[] | undefined = void 0;

    protected subscribe<T>(obs: Observable<T>, action: (v: T) => void) {
        if (typeof this.subs === 'undefined') this.subs = [];
        this.subs.push(obs.subscribe(action));
    }

    componentWillUnmount() {
        if (!this.subs) return;
        for (const s of this.subs) s.unsubscribe();
        this.subs = void 0;
    }

    protected init?(): void;

    constructor(props: P, context?: any) {
        super(props);
        this.plugin = context;
        if (this.init) this.init();
    }
}

export abstract class PurePluginUIComponent<P = {}, S = {}, SS = {}> extends React.PureComponent<P, S, SS> {
    static contextType = PluginReactContext;
    readonly plugin: PluginUIContext;

    private subs: Subscription[] | undefined = void 0;

    protected subscribe<T>(obs: Observable<T>, action: (v: T) => void) {
        if (typeof this.subs === 'undefined') this.subs = [];
        this.subs.push(obs.subscribe(action));
    }

    componentWillUnmount() {
        if (!this.subs) return;
        for (const s of this.subs) s.unsubscribe();
        this.subs = void 0;
    }

    protected init?(): void;

    constructor(props: P, context?: any) {
        super(props, context);
        this.plugin = context;
        if (this.init) this.init();
    }
}

export type _Props<C extends React.Component> = C extends React.Component<infer P> ? P : never
export type _State<C extends React.Component> = C extends React.Component<any, infer S> ? S : never
