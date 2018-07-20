/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const DefaultPlaneProps = {
    width: 1,
    height: 1
}
export type PlaneProps = Partial<typeof DefaultPlaneProps>

export function Plane(props?: PlaneProps) {
    const { width, height } = { ...DefaultPlaneProps, ...props }

    return {
        vertices: new Float32Array([
            -width / 2, height / 2, 0,
            width / 2, height / 2, 0,
            -width / 2, -height / 2, 0,
            width / 2, -height / 2, 0
        ]),
        normals: new Float32Array([
            0, 0, 1,
            0, 0, 1,
            0, 0, 1,
            0, 0, 1
        ]),
        indices: new Uint32Array([
            0, 2, 1,
            1, 2, 3
        ])
    }
}