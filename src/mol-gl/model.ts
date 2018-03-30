/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import { Mat4, Quat, Vec3 } from 'mol-math/linear-algebra'
import { defaults } from 'mol-util';

const tmpMat4 = Mat4()

type ModelProps = {
    rotation?: Quat,
    position?: Vec3,
    scale?: Vec3
}

function createModel(regl: REGL.Regl, props: ModelProps = {}) {
    const transform = Mat4.identity()

    const rotation = defaults(props.rotation, Quat.identity())
    const position = defaults(props.position, Vec3.zero())
    const scale = defaults(props.scale, Vec3.create(1, 1, 1))

    const draw = regl({
        context: { transform, rotation, position, scale },

        uniforms: {
            model(ctx: REGL.DefaultContext, props: any = {}) {
                const model = Mat4.identity()

                if ('rotation' in props) Quat.copy(rotation, props.rotation)
                if ('position' in props) Vec3.copy(position, props.position)
                if ('scale' in props) Vec3.copy(scale, props.scale)

                Mat4.translate(model, model, position)
                Mat4.mul(model, model, Mat4.fromQuat(tmpMat4, rotation))
                Mat4.scale(model, model, scale)

                if ('transform' in props) Mat4.mul(model, props.transform, model)
                Mat4.copy(transform, model)

                return model
            }
        }
    })

    return Object.assign(draw, {
        get transform() { return transform },
        get position() { return position },
        get rotation() { return rotation },
        get scale() { return scale },
    })
}

export default createModel