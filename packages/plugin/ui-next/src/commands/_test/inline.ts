/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { FieldHint, FormOption } from '../../forms/flat-fields.js';

/** A field in test shorthand: a text param, a string select (`options`), or a bounded number (`slider`). */
export interface TestField {
  id: string;
  label: string;
  before?: string;
  format?: string;
  options?: readonly FormOption[];
  optional?: boolean;
  slider?: { min: number; max: number; step?: number };
}

export function paramsOf(fields: readonly TestField[]): PD.Params {
  const params: PD.Params = {};
  for (const f of fields) {
    const param: PD.Any = f.options
      ? PD.Select(
          f.options[0]?.value as string,
          f.options.map((o) => [o.value, o.label] as const),
          { label: f.label },
        )
      : f.slider
        ? PD.Numeric(f.slider.min, f.slider, { label: f.label })
        : PD.Text('', { label: f.label });
    params[f.id] = f.optional ? (PD.Optional(param) as PD.Any) : param;
  }
  return params;
}

export function hintsOf(fields: readonly TestField[]): Record<string, FieldHint> {
  const hints: Record<string, FieldHint> = {};
  for (const f of fields) {
    if (f.before || f.format)
      hints[f.id] = { ...(f.before ? { before: f.before } : {}), ...(f.format ? { format: f.format } : {}) };
  }
  return hints;
}

/** The command parts for fixed inline fields: `...inline([...], 'Apply')`. */
export function inline(fields: readonly TestField[], action: string) {
  return { params: () => paramsOf(fields), fields: hintsOf(fields), action };
}
