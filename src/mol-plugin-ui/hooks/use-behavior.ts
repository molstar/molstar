/**
 * Copyright (c) 2020-21 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { useEffect, useRef, useState } from 'react';

interface Behavior<T> {
    value: T;
    subscribe(f: (v: T) => void): { unsubscribe(): void };
}

export function useBehavior<T>(s: Behavior<T>): T;
// eslint-disable-next-line
export function useBehavior<T>(s: Behavior<T> | undefined): T | undefined;
// eslint-disable-next-line
export function useBehavior<T>(s: Behavior<T> | undefined): T | undefined {
    const [, next] = useState({});
    const current = useRef<T>();
    current.current = s?.value;

    useEffect(() => {
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