/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginCommands } from 'mol-plugin/command';
import { PluginContext } from 'mol-plugin/context';
import { State, Transform } from 'mol-state';
import { StateAction } from 'mol-state/action';
import { memoizeOne } from 'mol-util/memoize';
import { StateTransformParameters, TransformContolBase } from './common';

export { ApplyActionContol };

namespace ApplyActionContol {
    export interface Props {
        plugin: PluginContext,
        nodeRef: Transform.Ref,
        state: State,
        action: StateAction
    }

    export interface ComponentState {
        ref: Transform.Ref,
        version: string,
        params: any,
        error?: string,
        busy: boolean,
        isInitial: boolean
    }
}

class ApplyActionContol extends TransformContolBase<ApplyActionContol.Props, ApplyActionContol.ComponentState> {
    applyAction() {
        return PluginCommands.State.ApplyAction.dispatch(this.plugin, {
            state: this.props.state,
            action: this.props.action.create(this.state.params),
            ref: this.props.nodeRef
        });
    }
    getInfo() { return this._getInfo(this.props.nodeRef, this.props.state.transforms.get(this.props.nodeRef).version); }
    getHeader() { return this.props.action.definition.display; }
    getHeaderFallback() { return this.props.action.id; }
    isBusy() { return !!this.state.error || this.state.busy; }
    applyText() { return 'Apply'; }

    private _getInfo = memoizeOne((t: Transform.Ref, v: string) => StateTransformParameters.infoFromAction(this.plugin, this.props.state, this.props.action, this.props.nodeRef));

    state = { ref: this.props.nodeRef, version: this.props.state.transforms.get(this.props.nodeRef).version, error: void 0, isInitial: true, params: this.getInfo().initialValues, busy: false };

    static getDerivedStateFromProps(props: ApplyActionContol.Props, state: ApplyActionContol.ComponentState) {
        if (props.nodeRef === state.ref) return null;
        const version = props.state.transforms.get(props.nodeRef).version;
        if (version === state.version) return null;

        const source = props.state.cells.get(props.nodeRef)!.obj!;
        const definition = props.action.definition.params || { };
        const initialValues = definition.default ? definition.default(source, props.plugin) : {};

        const newState: Partial<ApplyActionContol.ComponentState> = {
            ref: props.nodeRef,
            version,
            params: initialValues,
            isInitial: true,
            error: void 0
        };
        return newState;
    }
}