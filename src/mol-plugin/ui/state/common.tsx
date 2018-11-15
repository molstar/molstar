/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObject, State, Transform, StateObjectCell, Transformer } from 'mol-state';
import * as React from 'react';
import { PurePluginComponent } from '../base';
import { ParameterControls, ParamOnChange } from '../controls/parameters';
import { StateAction } from 'mol-state/action';
import { PluginContext } from 'mol-plugin/context';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Subject } from 'rxjs';

export { StateTransformParameters, TransformContolBase };

class StateTransformParameters extends PurePluginComponent<StateTransformParameters.Props> {
    validate(params: any) {
        // TODO
        return void 0;
    }

    areInitial(params: any) {
        return PD.areEqual(this.props.info.params, params, this.props.info.initialValues);
    }

    onChange: ParamOnChange = ({ name, value }) => {
        const params = { ...this.props.params, [name]: value };
        this.props.events.onChange(params, this.areInitial(params), this.validate(params));
    };

    render() {
        return <ParameterControls params={this.props.info.params} values={this.props.params} onChange={this.onChange} onEnter={this.props.events.onEnter} isDisabled={this.props.isDisabled} />;
    }
}


namespace StateTransformParameters {
    export interface Props {
        info: {
            params: PD.Params,
            initialValues: any,
            source: StateObject,
            isEmpty: boolean
        },
        events: {
            onChange: (params: any, areInitial: boolean, errors?: string[]) => void,
            onEnter: () => void,
        }
        params: any,
        isDisabled?: boolean
    }

    export type Class = React.ComponentClass<Props>

    export function infoFromAction(plugin: PluginContext, state: State, action: StateAction, nodeRef: Transform.Ref): Props['info'] {
        const source = state.cells.get(nodeRef)!.obj!;
        const params = action.definition.params ? action.definition.params(source, plugin) : { };
        const initialValues = PD.getDefaultValues(params);
        return {
            source,
            initialValues,
            params,
            isEmpty: Object.keys(params).length === 0
        };
    }

    export function infoFromTransform(plugin: PluginContext, state: State, transform: Transform): Props['info'] {
        const cell = state.cells.get(transform.ref)!;
        const source: StateObjectCell | undefined = (cell.sourceRef && state.cells.get(cell.sourceRef)!) || void 0;
        const create = transform.transformer.definition.params;
        const params = create ? create((source && source.obj) as any, plugin) : { };
        return {
            source: (source && source.obj) as any,
            initialValues: transform.params,
            params,
            isEmpty: Object.keys(params).length === 0
        }
    }
}

namespace TransformContolBase {
    export interface State {
        params: any,
        error?: string,
        busy: boolean,
        isInitial: boolean
    }
}

abstract class TransformContolBase<P, S extends TransformContolBase.State> extends PurePluginComponent<P, S> {
    abstract applyAction(): Promise<void>;
    abstract getInfo(): StateTransformParameters.Props['info'];
    abstract getHeader(): Transformer.Definition['display'];
    abstract getHeaderFallback(): string;
    abstract isBusy(): boolean;
    abstract applyText(): string;
    abstract state: S;

    private busy: Subject<boolean>;

    private onEnter = () => {
        if (this.state.error) return;
        this.apply();
    }

    events: StateTransformParameters.Props['events'] = {
        onEnter: this.onEnter,
        onChange: (params, isInitial, errors) => this.setState({ params, isInitial, error: errors && errors[0] })
    }

    apply = async () => {
        this.setState({ busy: true });
        try {
            await this.applyAction();
        } finally {
            this.busy.next(false);
        }
    }

    init() {
        this.busy = new Subject();
        this.subscribe(this.busy, busy => this.setState({ busy }));
    }

    refresh = () => {
        this.setState({ params: this.getInfo().initialValues, isInitial: true, error: void 0 });
    }

    setDefault = () => {
        const info = this.getInfo();
        const params = PD.getDefaultValues(info.params);
        this.setState({ params, isInitial: PD.areEqual(info.params, params, info.initialValues), error: void 0 });
    }

    render() {
        const info = this.getInfo();
        if (info.isEmpty) return <div>Nothing to update</div>;

        const display = this.getHeader();

        return <div>
            <div style={{ borderBottom: '1px solid #999', marginBottom: '5px' }}>
                <button onClick={this.setDefault} disabled={this.state.busy} style={{ float: 'right'}} title='Set default params'>↻</button>
                <h3>{(display && display.name) || this.getHeaderFallback()}</h3>
            </div>

            <StateTransformParameters info={info} events={this.events} params={this.state.params} isDisabled={this.state.busy} />

            <div style={{ textAlign: 'right' }}>
                <span style={{ color: 'red' }}>{this.state.error}</span>
                {this.state.isInitial ? void 0 : <button title='Refresh params' onClick={this.refresh} disabled={this.state.busy}>↶</button>}
                <button onClick={this.apply} disabled={this.isBusy()}>{this.applyText()}</button>
            </div>
        </div>
    }
}