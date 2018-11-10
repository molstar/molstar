/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginContext } from '../context';
import { Transform, State } from 'mol-state';
import { ParametersComponent } from 'mol-app/component/parameters';
import { StateAction } from 'mol-state/action';
import { PluginCommands } from 'mol-plugin/command';
import { UpdateTrajectory } from 'mol-plugin/state/actions/basic';

export class Controls extends React.Component<{ plugin: PluginContext }, { id: string }> {
    state = { id: '1grm' };

    private _snap: any = void 0;
    private getSnapshot = () => {
        this._snap = this.props.plugin.state.getSnapshot();
        console.log(btoa(JSON.stringify(this._snap)));
    }
    private setSnapshot = () => {
        if (!this._snap) return;
        this.props.plugin.state.setSnapshot(this._snap);
    }

    render() {
        return <div>
            <button onClick={() => this.props.plugin._test_centerView()}>Center View</button><br />
            <hr />
            <button onClick={this.getSnapshot}>Get Snapshot</button>
            <button onClick={this.setSnapshot}>Set Snapshot</button>
        </div>;
    }
}


export class _test_TrajectoryControls extends React.Component<{ plugin: PluginContext }> {
    render() {
        return <div>
            <b>Trajectory: </b>
            <button onClick={() => PluginCommands.State.ApplyAction.dispatch(this.props.plugin, {
                state: this.props.plugin.state.data,
                action: UpdateTrajectory.create({ action: 'advance', by: -1 })
            })}>&lt;&lt;</button>
            <button onClick={() => PluginCommands.State.ApplyAction.dispatch(this.props.plugin, {
                state: this.props.plugin.state.data,
                action: UpdateTrajectory.create({ action: 'reset' })
            })}>Reset</button>
            <button onClick={() => PluginCommands.State.ApplyAction.dispatch(this.props.plugin, {
                state: this.props.plugin.state.data,
                action: UpdateTrajectory.create({ action: 'advance', by: +1 })
            })}>&gt;&gt;</button><br />
        </div>
    }
}

export class _test_ApplyAction extends React.Component<{ plugin: PluginContext, nodeRef: Transform.Ref, state: State, action: StateAction }, { params: any }> {
    private getObj() {
        const obj = this.props.state.cells.get(this.props.nodeRef)!;
        return obj;
    }

    private getDefaultParams() {
        const p = this.props.action.definition.params;
        if (!p || !p.default) return {};
        const obj = this.getObj();
        if (!obj.obj) return {};
        return p.default(obj.obj, this.props.plugin);
    }

    private getParamDef() {
        const p = this.props.action.definition.params;
        if (!p || !p.controls) return {};
        const obj = this.getObj();
        if (!obj.obj) return {};
        return p.controls(obj.obj, this.props.plugin);
    }

    private create() {
        console.log('Apply Action', this.state.params);
        PluginCommands.State.ApplyAction.dispatch(this.props.plugin, {
            state: this.props.state,
            action: this.props.action.create(this.state.params),
            ref: this.props.nodeRef
        });
        // this.props.plugin.applyTransform(this.props.state, this.props.nodeRef, this.props.transformer, this.state.params);
    }

    state = { params: this.getDefaultParams() }

    render() {
        const obj = this.getObj();
        if (obj.status !== 'ok') {
            // TODO filter this elsewhere
            return <div />;
        }

        const action = this.props.action;

        return <div key={`${this.props.nodeRef} ${this.props.action.id}`}>
            <div style={{ borderBottom: '1px solid #999' }}><h3>{(action.definition.display && action.definition.display.name) || action.id}</h3></div>
            <ParametersComponent params={this.getParamDef()} values={this.state.params as any} onChange={(k, v) => {
                this.setState({ params: { ...this.state.params, [k]: v } });
            }} />
            <button onClick={() => this.create()} style={{ width: '100%' }}>Create</button>
        </div>
    }
}

export class _test_UpdateTransform extends React.Component<{ plugin: PluginContext, state: State, nodeRef: Transform.Ref }, { params: any }> {
    private getCell(ref?: string) {
        return this.props.state.cells.get(ref || this.props.nodeRef)!;
    }

    private getDefParams() {
        const cell = this.getCell();
        if (!cell) return {};
        return cell.transform.params;
    }

    private getParamDef() {
        const cell = this.getCell();
        const def = cell.transform.transformer.definition;

        if (!cell.sourceRef || !def.params || !def.params.controls) return void 0;
        const src = this.getCell(cell.sourceRef);
        if (!src || !src.obj) return void 0;

        return def.params.controls(src.obj, this.props.plugin);
    }

    private update() {
        console.log(this.props.nodeRef, this.state.params);
        this.props.plugin.updateTransform(this.props.state, this.props.nodeRef, this.state.params);
    }

    // componentDidMount() {
    //     const t = this.props.plugin.state.data.tree.nodes.get(this.props.nodeRef)!;
    //     if (t) this.setState({ params: t.value.params });
    // }

    state = { params: this.getDefParams() };

    render() {
        const cell = this.getCell();
        const transform = cell.transform;
        if (!transform || transform.ref === Transform.RootRef) {
            return <div />;
        }

        const params = this.getParamDef();
        if (!params) return <div />;

        const tr = transform.transformer;

        return <div key={`${this.props.nodeRef} ${tr.id}`} style={{ marginBottom: '10ox' }}>
            <div style={{ borderBottom: '1px solid #999' }}><h3>{(tr.definition.display && tr.definition.display.name) || tr.definition.name}</h3></div>
            <ParametersComponent params={params} values={this.state.params as any} onChange={(k, v) => {
                this.setState({ params: { ...this.state.params, [k]: v } });
            }} />
            <button onClick={() => this.update()} style={{ width: '100%' }}>Update</button>
        </div>
    }
}