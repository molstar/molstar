/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * The param definitions of value palettes (`value-palette.ts`): their values are `ValuePaletteProps`. The shape is
 * MolViewSpec's annotation color theme's (`mvs-annotation`, `palette`), which saved sessions hold: keep it.
 */

import { ParamDefinition as PD } from '../param-definition.js';
import { ColorNames } from './names.js';
import {
  type OptionalSplitColorProp as OptionalSplitColorPropT,
  type SplitColorProp as SplitColorPropT,
  SplitColorValue,
} from './split-color.js';
import type {
  CategoricalValuePaletteProps,
  ContinuousValuePaletteProps,
  DiscreteValuePaletteProps,
  ValuePaletteProps,
} from './value-palette.js';

/** Similar to `PD.Numeric` but allows leaving empty field in UI (treated as `undefined`) */
export function MaybeFloatParamDefinition(
  info?: PD.Info & { placeholder?: string },
): PD.Converted<number | null, string> {
  const defaultValue = null; // the default must be null, otherwise real nulls would be replaced by the default
  return PD.Converted(
    stringifyMaybeFloat,
    parseMaybeFloat,
    PD.Text(stringifyMaybeFloat(defaultValue), {
      ...info,
      disableInteractiveUpdates: true,
      placeholder: info?.placeholder ?? 'null',
    }),
  );
}
function parseMaybeFloat(input: string): number | null {
  const num = parseFloat(input);
  return isNaN(num) ? null : num;
}
function stringifyMaybeFloat(num: number | null): string {
  if (num === null) return '';
  return num.toString();
}

/** Return parameter definition for `SplitColorProp` */
function PDSplitColor(info?: PD.Info): PD.Mapped<SplitColorPropT> {
  return PD.MappedStatic(
    'oneColor',
    {
      oneColor: PD.Group(
        {
          color: PD.Color(ColorNames.white, { label: ' ' }),
        },
        { isFlat: true },
      ),
      splitColor: PD.Group(
        {
          color1: PD.Color(ColorNames.white, { label: '1' }),
          color2: PD.Color(ColorNames.white, { label: '2' }),
        },
        { isFlat: true },
      ),
    },
    info,
  );
}
/** Return parameter definition for `OptionalSplitColorProp` */
function PDOptionalSplitColor(info?: PD.Info): PD.Mapped<OptionalSplitColorPropT> {
  return PD.MappedStatic(
    'none',
    {
      none: PD.EmptyGroup(),
      oneColor: PD.Group(
        {
          color: PD.Color(ColorNames.white, { label: ' ' }),
        },
        { isFlat: true },
      ),
      splitColor: PD.Group(
        {
          color1: PD.Color(ColorNames.white, { label: '1' }),
          color2: PD.Color(ColorNames.white, { label: '2' }),
        },
        { isFlat: true },
      ),
    },
    info,
  );
}

export type SplitColorProp = SplitColorPropT;
export type OptionalSplitColorProp = OptionalSplitColorPropT;

/** Color parameter which allows regular color or split color (or optionally no color) */
export const SplitColorProp = {
  PD: PDSplitColor,
  PDOptional: PDOptionalSplitColor,
  ...SplitColorValue,
};

// The param objects' types are written out, so their values are the props `value-palette.ts` evaluates (inferred,
// the declaration files elide the object lists' element types).

export const CategoricalValuePaletteParams: {
  colors: PD.Mapped<CategoricalValuePaletteProps['colors']>;
  repeatColorList: PD.BooleanParam;
  sort: PD.Select<CategoricalValuePaletteProps['sort']>;
  sortDirection: PD.Select<CategoricalValuePaletteProps['sortDirection']>;
  caseInsensitive: PD.BooleanParam;
  missingColor: PD.Mapped<OptionalSplitColorPropT>;
} = {
  colors: PD.MappedStatic('list', {
    list: PD.ObjectList(
      {
        color: SplitColorProp.PD(),
      },
      (e) => SplitColorProp.toString(e.color),
      { description: 'List of colors.' },
    ),
    dictionary: PD.ObjectList(
      {
        value: PD.Text(),
        color: SplitColorProp.PD(),
      },
      (e) => `${e.value}: ${SplitColorProp.toString(e.color)}`,
      { description: 'Mapping of annotation values to colors.' },
    ),
  }),
  repeatColorList: PD.Boolean(false, {
    hideIf: (g) => g.colors.name !== 'list',
    description: 'Repeat color list once all colors are depleted (only applies if `colors` is a list).',
  }),
  sort: PD.Select(
    'none',
    [
      ['none', 'None'],
      ['lexical', 'Lexical'],
      ['numeric', 'Numeric'],
    ] as const,
    {
      hideIf: (g) => g.colors.name !== 'list',
      description:
        'Sort actual annotation values before assigning colors from a list (none = take values in order of their first occurrence).',
    },
  ),
  sortDirection: PD.Select(
    'ascending',
    [
      ['ascending', 'Ascending'],
      ['descending', 'Descending'],
    ] as const,
    { hideIf: (g) => g.colors.name !== 'list', description: 'Sort direction.' },
  ),
  caseInsensitive: PD.Boolean(false, { description: 'Treat annotation values as case-insensitive strings.' }),
  missingColor: SplitColorProp.PDOptional({
    description:
      'Color to use when (a) `colors` is a dictionary and given key is not present, or (b) `colors` is a list and there are more actual annotation values than listed colors and `repeat_color_list` is not true.',
  }),
};
export type CategoricalValuePaletteParams = typeof CategoricalValuePaletteParams;

export const DiscreteValuePaletteParams: {
  colors: PD.ObjectList<DiscreteValuePaletteProps['colors'][number]>;
  mode: PD.Select<DiscreteValuePaletteProps['mode']>;
  xMin: PD.Converted<number | null, string>;
  xMax: PD.Converted<number | null, string>;
} = {
  colors: PD.ObjectList(
    {
      color: SplitColorProp.PDOptional(),
      fromValue: PD.Numeric(-Infinity),
      toValue: PD.Numeric(Infinity),
    },
    (e) => `${SplitColorProp.toString(e.color)} [${formatFloat(e.fromValue)}, ${formatFloat(e.toValue)}]`,
    { description: 'Mapping of annotation value ranges to colors.' },
  ),
  mode: PD.Select(
    'normalized',
    [
      ['normalized', 'Normalized'],
      ['absolute', 'Absolute'],
    ] as const,
    {
      description:
        'Defines whether the annotation values should be normalized before assigning color based on checkpoints in `colors` (`x_normalized = (x - x_min) / (x_max - x_min)`, where `[x_min, x_max]` are either `value_domain` if provided, or the lowest and the highest value encountered in the annotation).',
    },
  ),
  xMin: MaybeFloatParamDefinition({
    hideIf: (g) => g.mode !== 'normalized',
    placeholder: 'auto',
    description:
      'Defines `x_min` for normalization of annotation values. If not provided, minimum of the actual values will be used. Only used when `mode` is `"normalized"',
  }),
  xMax: MaybeFloatParamDefinition({
    hideIf: (g) => g.mode !== 'normalized',
    placeholder: 'auto',
    description:
      'Defines `x_max` for normalization of annotation values. If not provided, maximum of the actual values will be used. Only used when `mode` is `"normalized"',
  }),
};
export type DiscreteValuePaletteParams = typeof DiscreteValuePaletteParams;

export const ContinuousValuePaletteParams: {
  colors: PD.ObjectList<ContinuousValuePaletteProps['colors'][number]>;
  mode: PD.Select<ContinuousValuePaletteProps['mode']>;
  xMin: PD.Converted<number | null, string>;
  xMax: PD.Converted<number | null, string>;
  underflowColor: PD.Mapped<OptionalSplitColorPropT>;
  overflowColor: PD.Mapped<OptionalSplitColorPropT>;
} = {
  colors: PD.ObjectList(
    {
      color: SplitColorProp.PD(),
      checkpoint: PD.Numeric(0),
    },
    (e) => `${SplitColorProp.toString(e.color)} [${formatFloat(e.checkpoint)}]`,
    { description: 'List of colors with checkpoints.' },
  ),
  mode: PD.Select(
    'normalized',
    [
      ['normalized', 'Normalized'],
      ['absolute', 'Absolute'],
    ] as const,
    {
      description:
        'Defines whether the annotation values should be normalized before assigning color based on checkpoints in `colors` (`x_normalized = (x - x_min) / (x_max - x_min)`, where `[x_min, x_max]` are either `value_domain` if provided, or the lowest and the highest value encountered in the annotation).',
    },
  ),
  xMin: MaybeFloatParamDefinition({
    hideIf: (g) => g.mode !== 'normalized',
    placeholder: 'auto',
    description:
      'Defines `x_min` for normalization of annotation values. If not provided, minimum of the actual values will be used. Only used when `mode` is `"normalized"',
  }),
  xMax: MaybeFloatParamDefinition({
    hideIf: (g) => g.mode !== 'normalized',
    placeholder: 'auto',
    description:
      'Defines `x_max` for normalization of annotation values. If not provided, maximum of the actual values will be used. Only used when `mode` is `"normalized"',
  }),
  underflowColor: SplitColorProp.PDOptional({ description: 'Color for values below the lowest checkpoint.' }),
  overflowColor: SplitColorProp.PDOptional({ description: 'Color for values above the highest checkpoint.' }),
};
export type ContinuousValuePaletteParams = typeof ContinuousValuePaletteParams;

/** A value palette: `direct` (the values are colors), categorical, discrete or continuous. */
export function ValuePaletteParam(info?: PD.Info): PD.Mapped<ValuePaletteProps> {
  return PD.MappedStatic(
    'direct',
    {
      direct: PD.EmptyGroup(),
      categorical: PD.Group(CategoricalValuePaletteParams),
      discrete: PD.Group(DiscreteValuePaletteParams),
      continuous: PD.Group(ContinuousValuePaletteParams),
    },
    info,
  );
}

function formatFloat(x: number): string {
  if (x === Infinity) return '∞';
  if (x === -Infinity) return '-∞';
  return x.toString();
}
