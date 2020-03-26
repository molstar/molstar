/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { Observable, Subscription } from 'rxjs';
import { PluginContext } from '../mol-plugin/context';
import { Icon } from './controls/icons';

export const PluginReactContext = React.createContext(void 0 as any as PluginContext);

export abstract class PluginUIComponent<P = {}, S = {}, SS = {}> extends React.Component<P, S, SS> {
    static contextType = PluginReactContext;
    readonly plugin: PluginContext;

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

export abstract class PurePluginUIComponent<P = {}, S = {}, SS = {}> extends React.PureComponent<P, S, SS> {
    static contextType = PluginReactContext;
    readonly plugin: PluginContext;

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

//

export type CollapsableProps = { initiallyCollapsed?: boolean, header?: string }
export type CollapsableState = { isCollapsed: boolean, header: string, description?: string, isHidden?: boolean }

export abstract class CollapsableControls<P = {}, S = {}, SS = {}> extends PluginUIComponent<P & CollapsableProps, S & CollapsableState, SS> {
    toggleCollapsed = () => {
        this.setState({ isCollapsed: !this.state.isCollapsed } as (S & CollapsableState))
    }

    protected abstract defaultState(): (S & CollapsableState)
    protected abstract renderControls(): JSX.Element | null

    render() {
        if (this.state.isHidden) return null;

        const wrapClass = this.state.isCollapsed
            ? 'msp-transform-wrapper msp-transform-wrapper-collapsed'
            : 'msp-transform-wrapper';

        return <div className={wrapClass}>
            <div className='msp-transform-header'>
                <button className='msp-btn msp-form-control msp-btn-block msp-btn-collapse msp-no-overflow' onClick={this.toggleCollapsed}>
                    <Icon name={this.state.isCollapsed ? 'expand' : 'collapse'} />
                    {this.state.header}
                    <small style={{ margin: '0 6px' }}>{this.state.isCollapsed ? '' : this.state.description}</small>
                </button>
            </div>
            {!this.state.isCollapsed && this.renderControls()}
        </div>
    }

    constructor(props: P & CollapsableProps, context?: any) {
        super(props, context)

        const state = this.defaultState()
        if (props.initiallyCollapsed !== undefined) state.isCollapsed = props.initiallyCollapsed
        if (props.header !== undefined) state.header = props.header
        this.state = state
    }
}