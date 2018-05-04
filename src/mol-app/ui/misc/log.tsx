/*
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol
 * Copyright (c) 2016 - now David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */

import * as React from 'react'
import { View } from '../view';
import { LogController } from '../../controller/misc/log';
import { CommonEvents } from '../../event/basic';
import { formatTime } from 'mol-util';
import { Logger } from '../../service/logger';

export class Log extends View<LogController, {}, {}> {

    private wrapper: HTMLDivElement | undefined = void 0;

    componentWillMount() {
        super.componentWillMount();
        this.subscribe(CommonEvents.LayoutChanged.getStream(this.controller.context), () => this.scrollToBottom());
    }

    componentDidUpdate() {
        this.scrollToBottom();
    }

    private scrollToBottom() {
        const log = this.wrapper;
        if (log) log.scrollTop = log.scrollHeight - log.clientHeight - 1;
    }

    render() {
        const entries = this.controller.latestState.entries;

        return <div className='molstar-log-wrap'>
            <div className='molstar-log' ref={log => this.wrapper = log!}>
                <ul className='molstar-list-unstyled'>
                    {entries.map((entry, i, arr) => {

                        let label: JSX.Element;
                        let e = entry!;
                        switch (e.type) {
                            case Logger.EntryType.Error:
                                label = <span className='label label-danger'>Error</span>;
                                break;
                            case Logger.EntryType.Warning:
                                label = <span className='label label-warning'>Warning</span>;
                                break;
                            case Logger.EntryType.Info:
                                label = <span className='label label-info'>Info</span>;
                                break;
                            default:
                                label = <span></span>
                        }

                        let t = formatTime(e.timestamp);
                        return <li key={i}>
                            <div className={'molstar-log-entry-badge molstar-log-entry-' + Logger.EntryType[e.type].toLowerCase()} />
                            {label}
                            <div className='molstar-log-timestamp'>{t}</div>
                            <div className='molstar-log-entry'>{e.message}</div>
                        </li>;
                    }) }
                </ul>
            </div>
        </div>;
    }
}