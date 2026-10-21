/**
 * Copyright (c) 2022-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { LinkBuilderProps, LinkCylinderProps } from '@molstar/graphics/repr/structure/visual/util/link';

export const InteractionsSharedParams = {
  sizeFactor: PD.Numeric(0.3, { min: 0, max: 10, step: 0.01 }),
  dashCount: PD.Numeric(6, { min: 1, max: 10, step: 1 }),
  dashScale: PD.Numeric(0.4, { min: 0, max: 2, step: 0.1 }),
  ignoreHydrogens: PD.Boolean(false),
  ignoreHydrogensVariant: PD.Select('all', PD.arrayToOptions(['all', 'non-polar'] as const)),
  tryUseImpostor: PD.Boolean(true),
  includeParent: PD.Boolean(false),
  parentDisplay: PD.Select('stub', PD.arrayToOptions(['stub', 'full', 'between'] as const), {
    description:
      'Only has an effect when "includeParent" is enabled. "Stub" shows just the child side of interactions to the parent. "Full" shows both sides of interactions to the parent. "Between" shows only interactions to the parent.',
  }),
};
export type InteractionsSharedParams = typeof InteractionsSharedParams;

/**
 * Props for `createLinkCylinderImpostors`: colorMode 'interpolate' needs a second location per group, which the
 * interactions location iterators do not give, so it falls back to 'default' (the mesh path ignores colorMode)
 */
export function impostorProps<P extends LinkCylinderProps>(props: P): P {
  return props.colorMode === 'default' ? props : { ...props, colorMode: 'default' };
}

/**
 * Dashed links with their impostor count, as many per link not ignored as `CylindersBuilder.addFixedCountDashes`
 * adds, so that
 * `createLinkCylinderImpostors` allocates its buffers once
 */
export function withDashedCylinderCount(builderProps: LinkBuilderProps, dashCount: number): LinkBuilderProps {
  const { linkCount, ignore } = builderProps;
  const cylinderCount = () => {
    let count = 0;
    for (let i = 0; i < linkCount; ++i) if (!ignore?.(i)) count++;
    return count * Math.floor((dashCount + 1) / 2);
  };
  return { ...builderProps, cylinderCount };
}
