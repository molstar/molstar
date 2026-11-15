/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec2 } from '@molstar/core/math/linear-algebra/3d/vec2';
import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ValueCell } from '@molstar/core/util';
import type { DataArray } from '@molstar/graphics/gl/renderable/util';
import { createDataArray } from '@molstar/graphics/gl/renderable/util';
import type { SiteColorValues } from '@molstar/graphics/gl/renderable/spheres';
import { ElementSymbolColors, elementSymbolColor } from '@molstar/graphics/theme/color/element-symbol';
import { type Structure, Unit } from '@molstar/model/model/structure';

// The second coloring of a particle target's atoms near the camera: by element, blended in
// by a sphere's radius on screen (`SiteColorSchema`). A sites target's by its sites' elements (`sites.ts`), a structure
// target's by its atoms' (here); its spheres' groups the elements, in the structure's serial order.

/** The second coloring's params, of the structure targets' spacefill (sites and structures) */
export const SiteColorParams = {
  siteColor: PD.MappedStatic(
    'element',
    {
      off: PD.EmptyGroup(),
      element: PD.Group({
        carbonColor: PD.Select<'particle' | 'element'>(
          'particle',
          [
            ['particle', 'Particle'],
            ['element', 'Element'],
          ],
          { description: "A carbon's color: its particle's (no shift of hue across the blend) or its element's." },
        ),
        radius: PD.Interval(
          [4, 12],
          { min: 0, max: 64, step: 0.5 },
          {
            description:
              'The radius of an atom on screen in pixels from which its element color blends in, and from which it is whole.',
          },
        ),
      }),
    },
    {
      description:
        'Atoms colored by element near the camera, blended in by their size on screen (atoms of known elements: a structure target, a particle pack of version 5, an mmCIF assembly).',
    },
  ),
};
export type SiteColorProps = PD.Values<typeof SiteColorParams>;
type SiteColorOn = Extract<SiteColorProps['siteColor'], { name: 'element' }>['params'];

/** The second color's values from its table: by the sites' types or by group, its range of radii on screen */
export function setSiteColorTable(
  values: SiteColorValues,
  table: DataArray<Uint8Array>,
  type: 'siteType' | 'group',
  params: SiteColorOn,
) {
  ValueCell.update(values.tSiteColor, table);
  ValueCell.updateIfChanged(values.dSiteColor, type);
  const [start, end] = params.radius;
  const range = values.uSiteColorRange.ref.value;
  if (range[0] !== start || range[1] !== end) ValueCell.update(values.uSiteColorRange, Vec2.set(range, start, end));
}

/**
 * The second color of a structure target's spheres (`SiteColorValues`): a texel an element (its group), its element's
 * color, a carbon's alpha 0 where it keeps the particle's color, as an element of no element (a coarse unit's) does;
 * off where asked
 */
export function setStructureTargetSiteColor(structure: Structure, props: SiteColorProps, values: SiteColorValues) {
  const p = props.siteColor;
  if (!p || p.name === 'off') {
    ValueCell.updateIfChanged(values.dSiteColor, 'off');
    return;
  }
  const keepCarbon = p.params.carbonColor === 'particle';
  const table = createDataArray(Math.max(1, structure.elementCount), 4, Uint8Array, values.tSiteColor.ref.value.array);
  const { array } = table;
  // in the structure's serial order, unit by unit, as its spheres' groups (`createStructureTargetSizes`)
  let g = 0;
  for (const unit of structure.units) {
    const { elements } = unit;
    if (!Unit.isAtomic(unit)) {
      array.fill(0, g * 4, (g + elements.length) * 4);
      g += elements.length;
      continue;
    }
    const symbols = unit.model.atomicHierarchy.atoms.type_symbol;
    for (let i = 0, il = elements.length; i < il; ++i, ++g) {
      const symbol = symbols.value(elements[i]);
      Color.toArray(elementSymbolColor(ElementSymbolColors, symbol), array, g * 4);
      array[g * 4 + 3] = keepCarbon && symbol === 'C' ? 0 : 255;
    }
  }
  setSiteColorTable(values, table, 'group', p.params);
}

/** The defines a structure target's impostors take of the second color (the precompile's signatures) */
export function getStructureSignatureDefines(props: SiteColorProps): { readonly dSiteColor: readonly string[] } {
  return { dSiteColor: [props.siteColor?.name === 'element' ? 'group' : 'off'] };
}
