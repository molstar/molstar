/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

const float maxFloat = 10000.0; // NOTE constant also set in in decodeFloatLog and in TypeScript
const float floatLogFactor = log(maxFloat + 1.0);
float encodeFloatLog(in float value) { return log(value + 1.0) / floatLogFactor; }

#pragma glslify: export(encodeFloatLog)