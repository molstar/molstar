/**
 * Copyright (c) 2020-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ColorNames } from '../../../mol-util/color/names';
import { Text } from '../../../mol-geo/geometry/text/text';

export const MeasurementRepresentationCommonTextParams = {
    customText: PD.Text('', { label: 'Text', description: 'Override the label with custom value.', isEssential: true }),
    textColor: PD.Color(ColorNames.black, { isEssential: true }),
    textSize: PD.Numeric(0.5, { min: 0.1, max: 10, step: 0.1 }, { isEssential: true }),
};

export const LociLabelTextParams = {
    ...Text.Params,
    ...MeasurementRepresentationCommonTextParams,
    borderWidth: PD.Numeric(0.2, { min: 0, max: 0.5, step: 0.01 })
};
export type LociLabelTextParams = typeof LociLabelTextParams