/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State, StateTransform } from 'mol-state';
import { memoizeLatest } from 'mol-util/memoize';
import { StateTransformParameters, TransformContolBase } from './common';
import { Observable } from 'rxjs';

export { UpdateTransformContol };

namespace UpdateTransformContol {
    export interface Props {
        transform: StateTransform,
        state: State,
        toggleCollapsed?: Observable<any>,
        initiallyCollapsed?: boolean
    }

    export interface ComponentState extends TransformContolBase.ComponentState {
        transform: StateTransform
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

    componentDidMount() {
        if (super.componentDidMount) super.componentDidMount();

        if (this.props.toggleCollapsed) this.subscribe(this.props.toggleCollapsed, () => this.setState({ isCollapsed: !this.state.isCollapsed }));

        this.subscribe(this.plugin.events.state.object.updated, ({ ref, state }) => {
            if (this.props.transform.ref !== ref || this.props.state !== state) return;
            if (this.state.params !== this.props.transform.params) {
                this._getInfo = memoizeLatest((t: StateTransform) => StateTransformParameters.infoFromTransform(this.plugin, this.props.state, t));
                this.setState({ params: this.props.transform.params, isInitial: true })
            }
        });
    }

    private _getInfo = memoizeLatest((t: StateTransform) => StateTransformParameters.infoFromTransform(this.plugin, this.props.state, t));

    state: UpdateTransformContol.ComponentState = { transform: this.props.transform, error: void 0, isInitial: true, params: this.getInfo().initialValues, busy: false, isCollapsed: this.props.initiallyCollapsed };

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