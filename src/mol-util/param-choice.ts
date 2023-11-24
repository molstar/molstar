/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { ParamDefinition } from './param-definition';


/**
 * Represents a set of values to choose from, with a default value. Example:
 * ```
 * export const MyChoice = new Choice({ yes: 'I agree', no: 'Nope' }, 'yes');
 * export type MyChoiceType = Choice.Values<typeof MyChoice>; // 'yes'|'no'
 * ```
 */
export class Choice<T extends string, D extends T> {
    readonly defaultValue: D;
    readonly options: [T, string][];
    private readonly nameDict: { [value in T]: string };
    constructor(opts: { [value in T]: string }, defaultValue: D) {
        this.defaultValue = defaultValue;
        this.options = Object.keys(opts).map(k => [k as T, opts[k as T]]);
        this.nameDict = opts;
    }
    PDSelect(defaultValue?: T, info?: ParamDefinition.Info): ParamDefinition.Select<T> {
        return ParamDefinition.Select<T>(defaultValue ?? this.defaultValue, this.options, info);
    }
    prettyName(value: T): string {
        return this.nameDict[value];
    }
    get values(): T[] {
        return this.options.map(([value, pretty]) => value);
    }
}
export namespace Choice {
    export type Values<T extends Choice<any, any>> = T extends Choice<infer R, any> ? R : any;
}
