/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State, StateTransform, StateTransformer } from '../../../mol-state';
import { memoizeLatest } from '../../../mol-util/memoize';
import { StateTransformParameters, TransformControlBase } from './common';
import { Observable } from 'rxjs';
import * as React from 'react';
import { PluginUIComponent } from '../base';

export { UpdateTransformControl, TransformUpdaterControl };

namespace UpdateTransformControl {
    export interface Props {
        transform: StateTransform,
        state: State,
        toggleCollapsed?: Observable<any>,
        initiallyCollapsed?: boolean,
        customHeader?: StateTransformer.Definition['display']
    }

    export interface ComponentState extends TransformControlBase.ComponentState {
        transform: StateTransform
    }
}

class UpdateTransformControl extends TransformControlBase<UpdateTransformControl.Props, UpdateTransformControl.ComponentState> {
    applyAction() { return this.plugin.updateTransform(this.props.state, this.props.transform.ref, this.state.params); }
    getInfo() { return this._getInfo(this.props.transform); }
    getTransformerId() { return this.props.transform.transformer.id; }
    getHeader() { return this.props.customHeader || this.props.transform.transformer.definition.display; }
    canApply() { return !this.state.error && !this.state.busy && !this.state.isInitial; }
    applyText() { return this.canApply() ? 'Update' : 'Nothing to Update'; }
    isUpdate() { return true; }
    getSourceAndTarget() {
        return {
            a: this.props.state.cells.get(this.props.transform.parent)!.obj,
            b: this.props.state.cells.has(this.props.transform.ref)! ? this.props.state.cells.get(this.props.transform.ref)!.obj : void 0
        };
    }

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

    state: UpdateTransformControl.ComponentState = { transform: this.props.transform, error: void 0, isInitial: true, params: this.getInfo().initialValues, busy: false, isCollapsed: this.props.initiallyCollapsed };

    static getDerivedStateFromProps(props: UpdateTransformControl.Props, state: UpdateTransformControl.ComponentState) {
        if (props.transform === state.transform) return null;
        const cell = props.state.cells.get(props.transform.ref)!;
        const newState: Partial<UpdateTransformControl.ComponentState> = {
            transform: props.transform,
            params: (cell.params && cell.params.values) || { },
            isInitial: true,
            error: void 0
        };
        return newState;
    }
}

class TransformUpdaterControl extends PluginUIComponent<{ nodeRef: string, initiallyCollapsed?: boolean, header?: StateTransformer.Definition['display'] }> {
    componentDidMount() {
        this.subscribe(this.plugin.events.state.object.updated, ({ ref, state }) => {
            if (this.props.nodeRef !== ref || this.plugin.state.dataState !== state) return;
            this.forceUpdate();
        });
    }

    render() {
        const state = this.plugin.state.dataState;
        const ref = this.props.nodeRef;
        const cell = state.cells.get(ref)!;

        if (!cell || (cell.status !== 'ok' && cell.status !== 'error')) return null;

        const transform = cell.transform;
        return <UpdateTransformControl state={state} transform={transform} initiallyCollapsed={this.props.initiallyCollapsed} customHeader={this.props.header} />;
    }
}