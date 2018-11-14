/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginCommands } from 'mol-plugin/command';
import { State, Transform } from 'mol-state';
import { StateAction } from 'mol-state/action';
import { Subject } from 'rxjs';
import { PurePluginComponent } from '../base';
import { StateTransformParameters } from './parameters';
import { memoizeOne } from 'mol-util/memoize';
import { PluginContext } from 'mol-plugin/context';

export { ApplyActionContol };

namespace ApplyActionContol {
    export interface Props {
        plugin: PluginContext,
        nodeRef: Transform.Ref,
        state: State,
        action: StateAction
    }

    export interface ComponentState {
        nodeRef: Transform.Ref,
        params: any,
        error?: string,
        busy: boolean,
        isInitial: boolean
    }
}

class ApplyActionContol extends PurePluginComponent<ApplyActionContol.Props, ApplyActionContol.ComponentState> {
    private busy: Subject<boolean>;

    onEnter = () => {
        if (this.state.error) return;
        this.apply();
    }

    source = this.props.state.cells.get(this.props.nodeRef)!.obj!;

    getInfo = memoizeOne((t: Transform.Ref) => StateTransformParameters.infoFromAction(this.plugin, this.props.state, this.props.action, this.props.nodeRef));

    events: StateTransformParameters.Props['events'] = {
        onEnter: this.onEnter,
        onChange: (params, isInitial, errors) => {
            this.setState({ params, isInitial, error: errors && errors[0] })
        }
    }

    // getInitialParams() {
    //     const p = this.props.action.definition.params;
    //     if (!p || !p.default) return {};
    //     return p.default(this.source, this.plugin);
    // }

    // initialErrors() {
    //     const p = this.props.action.definition.params;
    //     if (!p || !p.validate) return void 0;
    //     const errors = p.validate(this.info.initialValues, this.source, this.plugin);
    //     return errors && errors[0];
    // }

    state = { nodeRef: this.props.nodeRef, error: void 0, isInitial: true, params: this.getInfo(this.props.nodeRef).initialValues, busy: false };

    apply = async () => {
        this.setState({ busy: true });

        try {
            await PluginCommands.State.ApplyAction.dispatch(this.plugin, {
                state: this.props.state,
                action: this.props.action.create(this.state.params),
                ref: this.props.nodeRef
            });
        } finally {
            this.busy.next(false);
        }
    }

    init() {
        this.busy = new Subject();
        this.subscribe(this.busy, busy => this.setState({ busy }));
    }

    refresh = () => {
        this.setState({ params: this.getInfo(this.props.nodeRef).initialValues, isInitial: true, error: void 0 });
    }

    static getDerivedStateFromProps(props: ApplyActionContol.Props, state: ApplyActionContol.ComponentState) {
        if (props.nodeRef === state.nodeRef) return null;
        const source = props.state.cells.get(props.nodeRef)!.obj!;
        const definition = props.action.definition.params || { };
        const initialValues = definition.default ? definition.default(source, props.plugin) : {};

        const newState: Partial<ApplyActionContol.ComponentState> = {
            nodeRef: props.nodeRef,
            params: initialValues,
            isInitial: true,
            error: void 0
        };
        return newState;
    }

    render() {
        const info = this.getInfo(this.props.nodeRef);
        const action = this.props.action;

        return <div>
            <div style={{ borderBottom: '1px solid #999', marginBottom: '5px' }}><h3>{(action.definition.display && action.definition.display.name) || action.id}</h3></div>

            <StateTransformParameters info={info} events={this.events} params={this.state.params} isEnabled={!this.state.busy} />

            <div style={{ textAlign: 'right' }}>
                <span style={{ color: 'red' }}>{this.state.error}</span>
                {this.state.isInitial ? void 0 : <button title='Refresh Params' onClick={this.refresh} disabled={this.state.busy}>↻</button>}
                <button onClick={this.apply} disabled={!!this.state.error || this.state.busy}>Create</button>
            </div>
        </div>
    }
}