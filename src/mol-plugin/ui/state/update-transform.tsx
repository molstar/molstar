/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State, Transform } from 'mol-state';
import { memoizeLatest } from 'mol-util/memoize';
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
    getTransformerId() { return this.props.transform.transformer.id; }
    getHeader() { return this.props.transform.transformer.definition.display; }
    canApply() { return !this.state.error && !this.state.busy && !this.state.isInitial; }
    applyText() { return this.canApply() ? 'Update' : 'Nothing to Update'; }
    isUpdate() { return true; }

    canAutoApply(newParams: any) {
        const autoUpdate = this.props.transform.transformer.definition.canAutoUpdate
        if (!autoUpdate) return false;

        const { state } = this.props;
        const cell = state.cells.get(this.props.transform.ref);
        if (!cell || !cell.sourceRef || cell.status !== 'ok') return false;
        const parentCell = state.cells.get(cell.sourceRef)!;

        return autoUpdate({ a: cell.obj!, b: parentCell.obj!, oldParams: this.getInfo().initialValues, newParams }, this.plugin);
    }

    private _getInfo = memoizeLatest((t: Transform) => StateTransformParameters.infoFromTransform(this.plugin, this.props.state, this.props.transform));

    state: UpdateTransformContol.ComponentState = { transform: this.props.transform, error: void 0, isInitial: true, params: this.getInfo().initialValues, busy: false };

    static getDerivedStateFromProps(props: UpdateTransformContol.Props, state: UpdateTransformContol.ComponentState) {
        if (props.transform === state.transform) return null;
        const cell = props.state.cells.get(props.transform.ref)!;
        const newState: Partial<UpdateTransformContol.ComponentState> = {
            transform: props.transform,
            params: (cell.params && cell.params.values) || { },
            isInitial: true,
            error: void 0
        };
        return newState;
    }
}