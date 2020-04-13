/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/*
 * This code has been modified from https://github.com/toji/gl-matrix/,
 * copyright (c) 2015, Brandon Jones, Colin MacKenzie IV.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 */

import Mat4 from './3d/mat4';
import Mat3 from './3d/mat3';
import Vec2 from './3d/vec2';
import Vec3 from './3d/vec3';
import Vec4 from './3d/vec4';
import Quat from './3d/quat';
import { EPSILON } from './3d/common';

export { Mat4, Mat3, Vec2, Vec3, Vec4, Quat, EPSILON };

export type Vec<T> =
    T extends 4 ? Vec4 :
        T extends 3 ? Vec3 :
            T extends 2 ? Vec2 :
                number[]