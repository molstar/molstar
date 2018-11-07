/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginContext } from '../context';
import { Transform, Transformer, StateObject } from 'mol-state';
import { ParametersComponent } from 'mol-app/component/parameters';

export class Controls extends React.Component<{ plugin: PluginContext }, { id: string }> {
    state = { id: '1grm' };

    private createState = () => {
        const url = `http://www.ebi.ac.uk/pdbe/static/entry/${this.state.id.toLowerCase()}_updated.cif`;
        // const url = `https://webchem.ncbr.muni.cz/CoordinateServer/${this.state.id.toLowerCase()}/full`
        this.props.plugin._test_createState(url);
    }

    private _snap:any = void 0;
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
            <input type='text' defaultValue={this.state.id} onChange={e => this.setState({ id: e.currentTarget.value })} />
            <button onClick={this.createState}>Create State</button><br/>
            <button onClick={() => this.props.plugin._test_centerView()}>Center View</button><br/>
            <button onClick={() => this.props.plugin._test_nextModel()}>Next Model</button><br/>
            <button onClick={() => this.props.plugin._test_playModels()}>Play Models</button><br/>
            <hr />
            <button onClick={this.getSnapshot}>Get Snapshot</button>
            <button onClick={this.setSnapshot}>Set Snapshot</button>
        </div>;
    }
}

export class _test_CreateTransform extends React.Component<{ plugin: PluginContext, nodeRef: Transform.Ref, transformer: Transformer }, { params: any }> {
    private getObj() {
        const obj = this.props.plugin.state.data.objects.get(this.props.nodeRef)!;
        return obj;
    }

    private getDefaultParams() {
        const p = this.props.transformer.definition.params;
        if (!p || !p.default) return { };
        const obj = this.getObj();
        if (!obj.obj) return { };
        return p.default(obj.obj, this.props.plugin);
    }

    private getParamDef() {
        const p = this.props.transformer.definition.params;
        if (!p || !p.controls) return { };
        const obj = this.getObj();
        if (!obj.obj) return { };
        return p.controls(obj.obj, this.props.plugin);
    }

    private create() {
        console.log(this.props.transformer.definition.name, this.state.params);
        this.props.plugin._test_applyTransform(this.props.nodeRef, this.props.transformer, this.state.params);
    }

    state = { params: this.getDefaultParams() }

    render() {
        const obj = this.getObj();
        if (obj.state !== StateObject.StateType.Ok) {
            // TODO filter this elsewhere
            return <div />;
        }

        const t = this.props.transformer;

        return <div key={`${this.props.nodeRef} ${this.props.transformer.id}`}>
            <div style={{ borderBottom: '1px solid #999'}}>{(t.definition.display && t.definition.display.name) || t.definition.name}</div>
            <ParametersComponent params={this.getParamDef()} values={this.state.params as any} onChange={params => this.setState({ params })} />
            <button onClick={() => this.create()} style={{ width: '100%' }}>Create</button>
        </div>
    }
}