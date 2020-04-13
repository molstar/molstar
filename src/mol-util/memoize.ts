/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** Cache the latest result from calls to a function with any number of arguments */
export function memoizeLatest<Args extends any[], T>(f: (...args: Args) => T): (...args: Args) => T {
    let lastArgs: any[] | undefined = void 0, value: any = void 0;
    return (...args) => {
        if (!lastArgs || lastArgs.length !== args.length) {
            lastArgs = args;
            value = f.apply(void 0, args);
            return value;
        }
        for (let i = 0, _i = args.length; i < _i; i++) {
            if (args[i] !== lastArgs[i]) {
                lastArgs = args;
                value = f.apply(void 0, args);
                return value;
            }
        }
        return value;
    };
}

/** Cache all results from calls to a function with a single argument */
export function memoize1<A, T>(f: (a: A) => T): (a: A) => T {
    const cache = new Map<A, T>();
    return a => {
        if (cache.has(a)) return cache.get(a)!;
        const v = f(a);
        cache.set(a, v);
        return v;
    };
}