/**
 * Copyright (c) 2020-22 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import React from 'react';
import { skip } from 'rxjs';

interface Behavior<T> {
    value: T;
    subscribe(f: (v: T) => void): { unsubscribe(): void };
}

function useBehaviorLegacy<T>(s: Behavior<T> | undefined): T | undefined {
    const [, next] = React.useState({});
    const current = React.useRef<T>();
    current.current = s?.value;

    React.useEffect(() => {
        if (!s) {
            return;
        }
        const sub = s.subscribe((v) => {
            if (current.current !== v) next({});
        });

        return () => sub.unsubscribe();
    }, [s]);

    return s?.value;
}

function useBehaviorReact18<T>(s: Behavior<T> | undefined) {
    return (React as any).useSyncExternalStore(
        React.useCallback(
            (callback: () => void) => {
                const sub = (s as any)?.pipe!(skip(1)).subscribe(callback)!;
                return () => sub?.unsubscribe();
            },
            [s]
        ),
        React.useCallback(() => s?.value, [s])
    );
}

const _useBehavior = !!(React as any).useSyncExternalStore
    ? useBehaviorReact18
    : useBehaviorLegacy;

export function useBehavior<T>(s: Behavior<T>): T;
// eslint-disable-next-line
export function useBehavior<T>(s: Behavior<T> | undefined): T | undefined;
// eslint-disable-next-line
export function useBehavior<T>(s: Behavior<T> | undefined): T | undefined {
    return _useBehavior(s);
}