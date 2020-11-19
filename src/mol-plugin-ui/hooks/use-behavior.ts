/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { useEffect, useState } from 'react';

interface Behavior<T> {
    value: T;
    subscribe(f: (v: T) => void): { unsubscribe(): void };
}

export function useBehavior<T>(s: Behavior<T>): T;
// eslint-disable-next-line
export function useBehavior<T>(s: Behavior<T> | undefined): T | undefined;
// eslint-disable-next-line
export function useBehavior<T>(s: Behavior<T> | undefined): T | undefined {
    const [value, setValue] = useState(s?.value);

    useEffect(() => {
        if (!s) {
            if (value !== void 0) setValue(void 0);
            return;
        }
        let fst = true;
        const sub = s.subscribe((v) => {
            if (fst) {
                fst = false;
                if (v !== value) setValue(v);
            } else setValue(v);
        });

        return () => {
            sub.unsubscribe();
        };
        // eslint-disable-next-line
    }, [s]);

    return value;
}