/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State, Transform } from 'mol-state';
import * as React from 'react';
import { Subject } from 'rxjs';
import { PurePluginComponent } from '../base';
import { StateTransformParameters } from './parameters';
import { memoizeOne } from 'mol-util/memoize';

export { UpdateTransformContol };

namespace UpdateTransformContol {
    export interface Props {
        transform: Transform,
        state: State
    }

    export interface ComponentState {
        transform: Transform,
        params: any,
        error?: string,
        busy: boolean,
        isInitial: boolean
    }
}

class UpdateTransformContol extends PurePluginComponent<UpdateTransformContol.Props, UpdateTransformContol.ComponentState> {
    private busy: Subject<boolean>;

    onEnter = () => {
        if (this.state.error) return;
        this.apply();
    }

    getInfo = memoizeOne((t: Transform) => StateTransformParameters.infoFromTransform(this.plugin, this.props.state, this.props.transform));

    events: StateTransformParameters.Props['events'] = {
        onEnter: this.onEnter,
        onChange: (params, isInitial, errors) => {
            this.setState({ params, isInitial, error: errors && errors[0] })
        }
    }

    state: UpdateTransformContol.ComponentState = { transform: this.props.transform, error: void 0, isInitial: true, params: this.getInfo(this.props.transform).initialValues, busy: false };

    apply = async () => {
        this.setState({ busy: true });

        try {
            await this.plugin.updateTransform(this.props.state, this.props.transform.ref, this.state.params);
        } finally {
            this.busy.next(false);
        }
    }

    init() {
        this.busy = new Subject();
        this.subscribe(this.busy, busy => this.setState({ busy }));
    }

    refresh = () => {
        this.setState({ params: this.props.transform.params, isInitial: true, error: void 0 });
    }

    static getDerivedStateFromProps(props: UpdateTransformContol.Props, state: UpdateTransformContol.ComponentState) {
        if (props.transform === state.transform) return null;
        const newState: Partial<UpdateTransformContol.ComponentState> = {
            transform: props.transform,
            params: props.transform.params,
            isInitial: true,
            error: void 0
        };
        return newState;
    }

    render() {
        const info = this.getInfo(this.props.transform);
        if (info.isEmpty) return <div>Nothing to update</div>;

        const tr = this.props.transform.transformer;

        return <div>
            <div style={{ borderBottom: '1px solid #999', marginBottom: '5px' }}><h3>{(tr.definition.display && tr.definition.display.name) || tr.id}</h3></div>

            <StateTransformParameters info={info} events={this.events} params={this.state.params} isEnabled={!this.state.busy} />

            <div style={{ textAlign: 'right' }}>
                <span style={{ color: 'red' }}>{this.state.error}</span>
                {this.state.isInitial ? void 0 : <button title='Refresh Params' onClick={this.refresh} disabled={this.state.busy}>↻</button>}
                <button onClick={this.apply} disabled={!!this.state.error || this.state.busy || this.state.isInitial}>Update</button>
            </div>
        </div>
    }
}