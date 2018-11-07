/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginContext } from '../context';
import { StateTree } from './state-tree';
import { Viewport } from './viewport';
import { Controls, _test_CreateTransform } from './controls';
import { Transformer } from 'mol-state';

// TODO: base object with subscribe helpers

export class Plugin extends React.Component<{ plugin: PluginContext }, { }> {
    render() {
        return <div style={{ position: 'absolute', width: '100%', height: '100%' }}>
            <div style={{ position: 'absolute', width: '350px', height: '100%', overflowY: 'scroll' }}>
                <StateTree plugin={this.props.plugin} />
                <hr />
                <_test_CurrentObject plugin={this.props.plugin} />
            </div>
            <div style={{ position: 'absolute', left: '350px', right: '250px', height: '100%' }}>
                <Viewport plugin={this.props.plugin} />
            </div>
            <div style={{ position: 'absolute', width: '250px', right: '0', height: '100%' }}>
                <Controls plugin={this.props.plugin} />
            </div>
        </div>;
    }
}

export class _test_CurrentObject extends React.Component<{ plugin: PluginContext }, { }> {
    componentWillMount() {
        this.props.plugin.behaviors.state.data.currentObject.subscribe(() => this.forceUpdate());
    }
    render() {
        const ref = this.props.plugin.behaviors.state.data.currentObject.value.ref;
        // const n = this.props.plugin.state.data.tree.nodes.get(ref)!;
        const obj = this.props.plugin.state.data.objects.get(ref)!;

        const type = obj && obj.obj ? obj.obj.type : void 0;

        const transforms = type
            ? Transformer.fromType(type)
            : []
        return <div>
            Current Ref: {this.props.plugin.behaviors.state.data.currentObject.value.ref}
            <hr />
            {
                transforms.map((t, i) => <_test_CreateTransform key={`${t.id} ${ref} ${i}`} plugin={this.props.plugin} transformer={t} nodeRef={ref} />)
            }
        </div>;
    }
}