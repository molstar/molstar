/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObject, Transformer, State, Transform, StateObjectCell } from 'mol-state';
import { shallowEqual } from 'mol-util/object';
import * as React from 'react';
import { PurePluginComponent } from '../base';
import { ParameterControls, ParamOnChange } from '../controls/parameters';
import { StateAction } from 'mol-state/action';
import { PluginContext } from 'mol-plugin/context';
import { ParamDefinition as PD } from 'mol-util/param-definition';

export { StateTransformParameters };

class StateTransformParameters extends PurePluginComponent<StateTransformParameters.Props> {
    getDefinition() {
        const controls = this.props.info.definition.controls;
        if (!controls) return { };
        return controls!(this.props.info.source, this.plugin)
    }

    validate(params: any) {
        const validate = this.props.info.definition.validate;
        if (!validate) return void 0;
        return validate(params, this.props.info.source, this.plugin)
    }

    areInitial(params: any) {
        const areEqual = this.props.info.definition.areEqual;
        if (!areEqual) return shallowEqual(params, this.props.info.initialValues);
        return areEqual(params, this.props.info.initialValues);
    }

    onChange: ParamOnChange = ({ name, value }) => {
        const params = { ...this.props.params, [name]: value };
        this.props.events.onChange(params, this.areInitial(params), this.validate(params));
    };

    render() {
        return <ParameterControls params={this.props.info.params} values={this.props.params} onChange={this.onChange} onEnter={this.props.events.onEnter} isEnabled={this.props.isEnabled} />;
    }
}


namespace StateTransformParameters {
    export interface Props {
        info: {
            definition: Transformer.ParamsDefinition,
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
        isEnabled?: boolean
    }

    export type Class = React.ComponentClass<Props>

    export function infoFromAction(plugin: PluginContext, state: State, action: StateAction, nodeRef: Transform.Ref): Props['info'] {
        const source = state.cells.get(nodeRef)!.obj!;
        const definition = action.definition.params || { };
        const initialValues = definition.default ? definition.default(source, plugin) : {};
        const params = definition.controls ? definition.controls(source, plugin) : {};
        return {
            source,
            definition: action.definition.params || { },
            initialValues,
            params,
            isEmpty: Object.keys(params).length === 0
        };
    }

    export function infoFromTransform(plugin: PluginContext, state: State, transform: Transform): Props['info'] {
        const cell = state.cells.get(transform.ref)!;
        const source: StateObjectCell | undefined = (cell.sourceRef && state.cells.get(cell.sourceRef)!) || void 0;
        const definition = transform.transformer.definition.params || { };
        const params = definition.controls ? definition.controls((source && source.obj) as any, plugin) : {};
        return {
            source: (source && source.obj) as any,
            definition,
            initialValues: transform.params,
            params,
            isEmpty: Object.keys(params).length === 0
        }
    }
}