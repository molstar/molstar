/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

const float maxFloat = 10000.0; // NOTE constant also set in in encodeFloatLog and in TypeScript
const float floatLogFactor = log(maxFloat + 1.0);
float decodeFloatLog(in float value) { return exp(value * floatLogFactor) - 1.0; }

#pragma glslify: export(decodeFloatLog)