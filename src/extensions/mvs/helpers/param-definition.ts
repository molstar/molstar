/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';


/** Similar to `PD.Numeric` but allows leaving empty field in UI (treated as `undefined`) */
export function MaybeIntegerParamDefinition(defaultValue?: number, info?: PD.Info): PD.Base<number | undefined> {
    return PD.Converted<number | undefined, PD.Text>(stringifyMaybeInt, parseMaybeInt, PD.Text(stringifyMaybeInt(defaultValue), info));
}
/** The magic with negative zero looks crazy, but it's needed if we want to be able to write negative numbers, LOL. Please help if you know a better solution. */
function parseMaybeInt(input: string): number | undefined {
    if (input.trim() === '-') return -0;
    const num = parseInt(input);
    return isNaN(num) ? undefined : num;
}
function stringifyMaybeInt(num: number | undefined): string {
    if (num === undefined) return '';
    if (Object.is(num, -0)) return '-';
    return num.toString();
}

/** Similar to `PD.Text` but leaving empty field in UI is treated as `undefined` */
export function MaybeStringParamDefinition(defaultValue?: string, info?: PD.Info): PD.Base<string | undefined> {
    return PD.Converted<string | undefined, PD.Text>(stringifyMaybeString, parseMaybeString, PD.Text(stringifyMaybeString(defaultValue), info));
}
function parseMaybeString(input: string): string | undefined {
    return input === '' ? undefined : input;
}
function stringifyMaybeString(str: string | undefined): string {
    return str === undefined ? '' : str;
}
