/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export { LogEntry };

interface LogEntry {
    type: LogEntry.Type,
    timestamp: Date,
    message: string
}

namespace LogEntry {
    export type Type = 'message' | 'error' | 'warning' | 'info'

    export function message(msg: string): LogEntry { return { type: 'message', timestamp: new Date(), message: msg }; }
    export function error(msg: string): LogEntry { return { type: 'error', timestamp: new Date(), message: msg }; }
    export function warning(msg: string): LogEntry { return { type: 'warning', timestamp: new Date(), message: msg }; }
    export function info(msg: string): LogEntry { return { type: 'info', timestamp: new Date(), message: msg }; }
}