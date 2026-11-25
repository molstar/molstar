/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color, ColorScale } from '@molstar/core/util/color';
import type { Location } from '@molstar/model/model/location';
import type { ColorTheme, LocationColor } from '../color.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ThemeDataContext } from '../theme.js';
import type { ScaleLegend, TableLegend } from '@molstar/core/util/legend';
import { ColorThemeCategory } from './categories.js';
import { fillByLocation, getPagedPainting, setParticleRgba } from './particle-pages.js';
import { Particle, type ParticleList } from '@molstar/model/model/particles/particle-list';

const DefaultColor = Color(0xcccccc);
const Description = 'Colors particles by a per-particle scalar attribute (e.g. score, class, CC).';

function getParticleList(ctx: ThemeDataContext): ParticleList | undefined {
  if (ctx.particles) return ctx.particles;
  return undefined;
}

/** A list's attributes: a slot list's its first page's (a tile set's pages have the same, the whole's ranges) */
function getAttributes(particles: ParticleList | undefined) {
  return particles?.attributes ?? particles?.pages?.pages[0]?.list.attributes;
}

export const ParticleAttributeColorThemeParams = {
  attribute: PD.Select('', [['', 'None']] as [string, string][]),
  domain: PD.MappedStatic('auto', {
    auto: PD.Group({}),
    custom: PD.Interval([0, 1], { step: 0.001 }),
  }),
  list: PD.ColorList('red-white-blue', { presetKind: 'scale' }),
};
export type ParticleAttributeColorThemeParams = typeof ParticleAttributeColorThemeParams;

export function getParticleAttributeColorThemeParams(ctx: ThemeDataContext) {
  const params = PD.clone(ParticleAttributeColorThemeParams);
  const particles = getParticleList(ctx);
  const attributes = getAttributes(particles);
  if (attributes && attributes.size > 0) {
    const options: [string, string][] = [];
    attributes.forEach((attr, key) => {
      if (['int', 'float'].includes(attr.column.schema.valueType)) {
        options.push([key, attr.label]);
      }
    });
    if (options.length > 0) {
      params.attribute = PD.Select(options[0][0], options);
    }
  }
  return params;
}

export function ParticleAttributeColorTheme(
  ctx: ThemeDataContext,
  props: PD.Values<ParticleAttributeColorThemeParams>,
): ColorTheme<ParticleAttributeColorThemeParams> {
  let color: LocationColor;
  let legend: ScaleLegend | TableLegend | undefined;
  let fillParticleColors: ColorTheme<ParticleAttributeColorThemeParams>['fillParticleColors'];

  const particles = getParticleList(ctx);
  const attributes = getAttributes(particles);
  if (attributes) {
    const attr = attributes.get(props.attribute);

    if (attr) {
      const domain: [number, number] =
        props.domain.name === 'custom' ? (props.domain.params as [number, number]) : [attr.min, attr.max];

      const scale = ColorScale.create({
        reverse: false,
        domain,
        listOrName: props.list.colors,
      });
      legend = scale.legend;

      const pick = (index: number): Color => {
        if (index < 0 || index >= attr.column.rowCount) return DefaultColor;
        return scale.color(attr.column.value(index));
      };

      if (particles!.pages) {
        // a slot list: a page's own attribute
        const painting = getPagedPainting((page) => {
          const column = page.attributes?.get(props.attribute)?.column;
          const at = (i: number) =>
            !column || i < 0 || i >= column.rowCount ? DefaultColor : scale.color(column.value(i));
          return {
            color: at,
            fill: (out, offset) => {
              for (let i = 0; i < page.count; ++i) setParticleRgba(out, offset + i * 4, at(i));
            },
          };
        }, DefaultColor);
        color = painting.color;
        fillParticleColors = painting.fill;
      } else {
        color = (location: Location): Color => {
          if (Particle.isLocation(location)) {
            return pick(location.index);
          }
          return DefaultColor;
        };
        fillParticleColors = (list, out, offset) => {
          if (list !== particles) return fillByLocation(color, list, out, offset);
          for (let i = 0; i < list.count; ++i) setParticleRgba(out, offset + i * 4, pick(i));
        };
      }
    } else {
      color = () => DefaultColor;
    }
  } else {
    color = () => DefaultColor;
  }

  return {
    factory: ParticleAttributeColorTheme,
    granularity: 'instance',
    color,
    props,
    description: Description,
    legend,
    fillParticleColors,
  };
}

export const ParticleAttributeColorThemeProvider: ColorTheme.Provider<
  ParticleAttributeColorThemeParams,
  'particle-attribute'
> = {
  name: 'particle-attribute',
  label: 'Particle Attribute',
  category: ColorThemeCategory.Particle,
  factory: ParticleAttributeColorTheme,
  getParams: getParticleAttributeColorThemeParams,
  defaultValues: PD.getDefaultValues(ParticleAttributeColorThemeParams),
  isApplicable: (ctx: ThemeDataContext) => {
    const attributes = getAttributes(getParticleList(ctx));
    return !!attributes && attributes.size > 0;
  },
};
