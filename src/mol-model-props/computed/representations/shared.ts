/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';

export const InteractionsSharedParams = {
    sizeFactor: PD.Numeric(0.3, { min: 0, max: 10, step: 0.01 }),
    dashCount: PD.Numeric(6, { min: 1, max: 10, step: 1 }),
    dashScale: PD.Numeric(0.4, { min: 0, max: 2, step: 0.1 }),
    ignoreHydrogens: PD.Boolean(false),
    ignoreHydrogensVariant: PD.Select('all', PD.arrayToOptions(['all', 'non-polar'] as const)),
    includeParent: PD.Boolean(false),
    parentDisplay: PD.Select('stub', PD.arrayToOptions(['stub', 'full', 'between'] as const), { description: 'Only has an effect when "includeParent" is enabled. "Stub" shows just the child side of interactions to the parent. "Full" shows both sides of interactions to the parent. "Between" shows only interactions to the parent.' }),
};
export type InteractionsSharedParams = typeof InteractionsSharedParams
