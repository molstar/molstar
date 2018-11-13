/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { Transform, State, Transformer } from 'mol-state';
import { StateAction } from 'mol-state/action';
import { PluginCommands } from 'mol-plugin/command';
import { PluginComponent } from './base';
import { ParameterControls, createChangeSubject, ParamChanges } from './controls/parameters';
import { Subject } from 'rxjs';
import { shallowEqual } from 'mol-util/object';

export { ActionContol }

namespace ActionContol {
    export interface Props {
        nodeRef: Transform.Ref,
        state: State,
        action?: StateAction
    }
}

class ActionContol extends PluginComponent<ActionContol.Props, { params: any, initialParams: any, error?: string, busy: boolean, canApply: boolean }> {
    private changes: ParamChanges;
    private busy: Subject<boolean>;

    cell = this.props.state.cells.get(this.props.nodeRef)!;
    parentCell = (this.cell.sourceRef && this.props.state.cells.get(this.cell.sourceRef)) || void 0;

    action: StateAction | Transformer = !this.props.action ? this.cell.transform.transformer : this.props.action
    isUpdate = !this.props.action

    getDefaultParams() {
        if (this.isUpdate) {
            return this.cell.transform.params;
        } else {
            const p = this.action.definition.params;
            if (!p || !p.default) return {};
            const obj = this.cell;
            if (!obj.obj) return {};
            return p.default(obj.obj, this.plugin);
        }
    }

    getParamDefinitions() {
        if (this.isUpdate) {
            const cell = this.cell;
            const def = cell.transform.transformer.definition;

            if (!cell.sourceRef || !def.params || !def.params.controls) return { };
            const src = this.parentCell;
            if (!src || !src.obj) return { };

            return def.params.controls(src.obj, this.plugin);
        } else {
            const p = this.action.definition.params;
            if (!p || !p.controls) return {};
            const cell = this.cell;
            if (!cell.obj) return {};
            return p.controls(cell.obj, this.plugin);
        }
    }

    defaultState() {
        const params = this.getDefaultParams();
        return { error: void 0, params, initialParams: params, busy: false, canApply: !this.isUpdate };
    }

    apply = async () => {
        this.setState({ busy: true, initialParams: this.state.params, canApply: !this.isUpdate });

        try {
            if (Transformer.is(this.action)) {
                await this.plugin.updateTransform(this.props.state, this.props.nodeRef, this.state.params);
            } else {
                await PluginCommands.State.ApplyAction.dispatch(this.plugin, {
                    state: this.props.state,
                    action: this.action.create(this.state.params),
                    ref: this.props.nodeRef
                });
            }
        } finally {
            this.busy.next(false);
        }
    }

    validate(params: any) {
        const def = this.isUpdate ? this.cell.transform.transformer.definition.params : this.action.definition.params;
        if (!def || !def.validate) return;
        const cell = this.cell;
        const error = def.validate(params, this.isUpdate ? this.parentCell!.obj! : cell.obj!, this.plugin);
        return error && error[0];
    }

    init() {
        this.changes = createChangeSubject();
        this.subscribe(this.changes, ({ name, value }) => {
            const params = { ...this.state.params, [name]: value };
            const canApply = this.isUpdate ? !shallowEqual(params, this.state.initialParams) : true;
            this.setState({ params, error: this.validate(params), canApply });
        });

        this.busy = new Subject();
        this.subscribe(this.busy, busy => this.setState({ busy }));
    }

    onEnter = () => {
        if (this.state.error) return;
        this.apply();
    }

    refresh = () => {
        this.setState({ params: this.state.initialParams, canApply: !this.isUpdate });
    }

    state = this.defaultState()

    render() {
        console.log('render', this.props.nodeRef, this.action.id);
        const cell = this.cell;
        if (cell.status !== 'ok' || (this.isUpdate && cell.transform.ref === Transform.RootRef)) return null;

        const action = this.action;

        return <div>
            <div style={{ borderBottom: '1px solid #999', marginBottom: '5px' }}><h3>{(action.definition.display && action.definition.display.name) || action.id}</h3></div>

            <ParameterControls params={this.getParamDefinitions()} values={this.state.params} changes={this.changes} onEnter={this.onEnter} isEnabled={!this.state.busy} />

            <div style={{ textAlign: 'right' }}>
                <span style={{ color: 'red' }}>{this.state.error}</span>
                <button onClick={this.apply} disabled={!this.state.canApply || !!this.state.error || this.state.busy}>{this.isUpdate ? 'Update' : 'Create'}</button>
                <button title='Refresh Params' onClick={this.refresh} disabled={this.state.busy}>↻</button>
            </div>
        </div>
    }
}