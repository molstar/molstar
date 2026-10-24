/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

/** Similar to `PD.Numeric` but allows leaving empty field in UI (treated as `undefined`), and can only have integer values */
export function MaybeIntegerParamDefinition(
  info?: PD.Info & { placeholder?: string },
): PD.Converted<number | null, string> {
  const defaultValue = null; // the default must be null, otherwise real nulls would be replaced by the default
  return PD.Converted(
    stringifyMaybeInt,
    parseMaybeInt,
    PD.Text(stringifyMaybeInt(defaultValue), {
      ...info,
      disableInteractiveUpdates: true,
      placeholder: info?.placeholder ?? 'null',
    }),
  );
}
function parseMaybeInt(input: string): number | null {
  const num = parseInt(input);
  return isNaN(num) ? null : num;
}
function stringifyMaybeInt(num: number | null): string {
  if (num === null) return '';
  return num.toString();
}

export { MaybeFloatParamDefinition } from '@molstar/core/util/color/value-palette-params';

/** Similar to `PD.Text` but leaving empty field in UI is treated as `undefined` */
export function MaybeStringParamDefinition(
  info?: PD.Info & { placeholder?: string },
): PD.Converted<string | null, string> {
  const defaultValue = null; // the default must be null, otherwise real nulls would be replaced by the default
  return PD.Converted(
    stringifyMaybeString,
    parseMaybeString,
    PD.Text(stringifyMaybeString(defaultValue), { ...info, placeholder: info?.placeholder ?? 'null' }),
  );
}
function parseMaybeString(input: string): string | null {
  return input === '' ? null : input;
}
function stringifyMaybeString(str: string | null): string {
  return str === null ? '' : str;
}
