/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol (c) David Sehnal
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginUIComponent } from './base';
import { PluginToastManager } from '../mol-plugin/util/toast';
import { IconButton } from './controls/common';
import { CancelSvg } from './controls/icons';

class ToastEntry extends PluginUIComponent<{ entry: PluginToastManager.Entry }> {
    private hide = () => {
        let entry = this.props.entry;
        (entry.hide || function () { }).call(null);
    };

    render() {
        let entry = this.props.entry;
        let message = typeof entry.message === 'string'
            ? <div dangerouslySetInnerHTML={{ __html: entry.message }} />
            : <div><entry.message /></div>;

        return <div className='msp-toast-entry'>
            <div className='msp-toast-title' onClick={() => this.hide()}>
                {entry.title}
            </div>
            <div className='msp-toast-message'>
                {message}
            </div>
            <div className='msp-toast-clear'></div>
            <div className='msp-toast-hide'>
                <IconButton svg={CancelSvg} onClick={this.hide} title='Hide' className='msp-no-hover-outline' />
            </div>
        </div>;
    }
}

export class Toasts extends PluginUIComponent {
    componentDidMount() {
        this.subscribe(this.plugin.managers.toast.events.changed, () => this.forceUpdate());
    }

    render() {
        const state = this.plugin.managers.toast.state;

        if (!state.entries.count()) return null;

        const entries: PluginToastManager.Entry[] = [];
        state.entries.forEach((t, k) => entries.push(t!));
        entries.sort(function (x, y) { return x.serialNumber - y.serialNumber; });

        return <div className='msp-toast-container'>
            {entries.map(e => <ToastEntry key={e.serialNumber} entry={e} />)}
        </div>;
    }
}
