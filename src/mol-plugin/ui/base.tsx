/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { Observable, Subscription } from 'rxjs';
import { PluginContext } from '../context';

export const PluginReactContext = React.createContext(void 0 as any as PluginContext);

export abstract class PluginUIComponent<P = {}, S = {}, SS = {}> extends React.Component<P, S, SS> {
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

export abstract class PurePluginUIComponent<P = {}, S = {}, SS = {}> extends React.PureComponent<P, S, SS> {
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

export type _Props<C extends React.Component> = C extends React.Component<infer P> ? P : never
export type _State<C extends React.Component> = C extends React.Component<any, infer S> ? S : never

//

export type CollapsableProps = { initiallyCollapsed?: boolean, header?: string }
export type CollapsableState = { isCollapsed: boolean, header: string }

export abstract class CollapsableControls<P extends CollapsableProps = CollapsableProps, S extends CollapsableState = CollapsableState, SS = {}> extends PluginUIComponent<P, S, SS> {
    toggleCollapsed = () => {
        this.setState({ isCollapsed: !this.state.isCollapsed })
    }

    protected abstract defaultState(): S
    protected abstract renderControls(): JSX.Element | null

    render() {
        const wrapClass = this.state.isCollapsed
            ? 'msp-transform-wrapper msp-transform-wrapper-collapsed'
            : 'msp-transform-wrapper';

        return <div className={wrapClass}>
            <div className='msp-transform-header'>
                <button className='msp-btn msp-btn-block' onClick={this.toggleCollapsed}>
                    <span className={`msp-icon msp-icon-${this.state.isCollapsed ? 'expand' : 'collapse'}`} />
                    {this.state.header}
                </button>
            </div>
            {!this.state.isCollapsed && this.renderControls()}
        </div>
    }

    constructor(props: P, context?: any) {
        super(props, context)
        this.state = this.defaultState()
    }
}