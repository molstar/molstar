/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State, Transform } from 'mol-state';
import { memoizeOne } from 'mol-util/memoize';
import { StateTransformParameters, TransformContolBase } from './common';

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

class UpdateTransformContol extends TransformContolBase<UpdateTransformContol.Props, UpdateTransformContol.ComponentState> {
    applyAction() { return this.plugin.updateTransform(this.props.state, this.props.transform.ref, this.state.params); }
    getInfo() { return this._getInfo(this.props.transform); }
    getHeader() { return this.props.transform.transformer.definition.display; }
    getHeaderFallback() { return this.props.transform.transformer.definition.name; }
    isBusy() { return !!this.state.error || this.state.busy || this.state.isInitial; }
    applyText() { return 'Update'; }

    private _getInfo = memoizeOne((t: Transform) => StateTransformParameters.infoFromTransform(this.plugin, this.props.state, this.props.transform));

    state: UpdateTransformContol.ComponentState = { transform: this.props.transform, error: void 0, isInitial: true, params: this.getInfo().initialValues, busy: false };

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
}