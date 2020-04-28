/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ColorNames } from '../../../mol-util/color/names';

export const MeasurementRepresentationCommonTextParams = {
    customText: PD.Text('', { label: 'Text', description: 'Override the label with custom value.' }),
    textColor: PD.Color(ColorNames.black, { isEssential: true }),
    textSize: PD.Numeric(0.5, { min: 0.1, max: 5, step: 0.1 }, { isEssential: true }),
};