/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export function memoizeOne<Args extends any[], T>(f: (...args: Args) => T): (...args: Args) => T {
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
    }
}