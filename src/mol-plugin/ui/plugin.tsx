/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginContext } from '../context';
import { Tree } from './tree';
import { Viewport } from './viewport';
import { Controls } from './controls';

export class Plugin extends React.Component<{ plugin: PluginContext }, { }> {
    render() {
        return <div style={{ position: 'absolute', width: '100%', height: '100%' }}>
            <div style={{ position: 'absolute', width: '250px', height: '100%' }}>
                <Tree plugin={this.props.plugin} />
            </div>
            <div style={{ position: 'absolute', left: '250px', right: '250px', height: '100%' }}>
                <Viewport plugin={this.props.plugin} />
            </div>
            <div style={{ position: 'absolute', width: '250px', right: '0', height: '100%' }}>
                <Controls plugin={this.props.plugin} />
            </div>
        </div>;
    }
}