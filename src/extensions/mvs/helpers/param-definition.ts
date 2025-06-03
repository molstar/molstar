/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';


/** Similar to `PD.Numeric` but allows leaving empty field in UI (treated as `undefined`), and can only have integer values */
export function MaybeIntegerParamDefinition(defaultValue?: number, info?: PD.Info & { placeholder?: string }): PD.Converted<number | undefined, string> {
    return PD.Converted(stringifyMaybeInt, parseMaybeInt, PD.Text(stringifyMaybeInt(defaultValue), { ...info, disableInteractiveUpdates: true, placeholder: info?.placeholder ?? 'undefined' }));
}
function parseMaybeInt(input: string): number | undefined {
    const num = parseInt(input);
    return isNaN(num) ? undefined : num;
}
function stringifyMaybeInt(num: number | undefined): string {
    if (num === undefined) return '';
    return num.toString();
}


/** Similar to `PD.Numeric` but allows leaving empty field in UI (treated as `undefined`) */
export function MaybeFloatParamDefinition(defaultValue?: number, info?: PD.Info & { placeholder?: string }): PD.Converted<number | undefined, string> {
    return PD.Converted(stringifyMaybeFloat, parseMaybeFloat, PD.Text(stringifyMaybeFloat(defaultValue), { ...info, disableInteractiveUpdates: true, placeholder: info?.placeholder ?? 'undefined' }));
}
function parseMaybeFloat(input: string): number | undefined {
    const num = parseFloat(input);
    return isNaN(num) ? undefined : num;
}
function stringifyMaybeFloat(num: number | undefined): string {
    if (num === undefined) return '';
    return num.toString();
}


/** Similar to `PD.Text` but leaving empty field in UI is treated as `undefined` */
export function MaybeStringParamDefinition(defaultValue?: string, info?: PD.Info): PD.Converted<string | undefined, string> {
    return PD.Converted(stringifyMaybeString, parseMaybeString, PD.Text(stringifyMaybeString(defaultValue), info));
}
function parseMaybeString(input: string): string | undefined {
    return input === '' ? undefined : input;
}
function stringifyMaybeString(str: string | undefined): string {
    return str === undefined ? '' : str;
}
