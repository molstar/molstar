/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

// TODO workaround due to some kind of GPU quirk
float myDiv(float a, float b) {
    return float(int(a) / int(b));
}

#pragma glslify: export(myDiv)