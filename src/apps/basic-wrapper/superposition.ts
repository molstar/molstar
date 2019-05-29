/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// TODO: move to an "example"

import { PluginContext } from 'mol-plugin/context';
import { Mat4 } from 'mol-math/linear-algebra';
import { StateHelper } from './helpers';

export type SuperpositionTestInput = {
    pdbId: string,
    auth_asym_id: string,
    matrix: Mat4
}[];

// function getAxisAngleTranslation(m: Mat4) {
//     const translation = Mat4.getTranslation(Vec3.zero(), m);
//     const axis = Vec3.zero();
//     const angle = 180 / Math.PI * Quat.getAxisAngle(axis, Mat4.getRotation(Quat.zero(), m));
//     return { translation, axis, angle };
// }

export function buildStaticSuperposition(ctx: PluginContext, src: SuperpositionTestInput) {
    const b = ctx.state.dataState.build().toRoot();
    for (const s of src) {
        StateHelper.visual(ctx,
            StateHelper.transform(
                StateHelper.selectChain(
                    StateHelper.structure(
                        StateHelper.getModel(StateHelper.download(b, `https://www.ebi.ac.uk/pdbe/static/entry/${s.pdbId}_updated.cif`), 'cif')),
                    s.auth_asym_id
                ),
                s.matrix
            )
        );
    }
    return b;
}

export const StaticSuperpositionTestData: SuperpositionTestInput = [
    { pdbId: '1aj5', auth_asym_id: 'A', matrix: Mat4.identity() },
    { pdbId: '1df0', auth_asym_id: 'B', matrix: Mat4.ofRows(
        [[0.406, 0.879, 0.248, -200.633],
         [0.693, -0.473, 0.544, 73.403],
         [0.596, -0.049, -0.802, -14.209],
         [0, 0, 0, 1]] )},
    { pdbId: '1dvi', auth_asym_id: 'A', matrix: Mat4.ofRows(
        [[-0.053, -0.077, 0.996, -45.633],
         [-0.312, 0.949, 0.057, -12.255],
         [-0.949, -0.307, -0.074, 53.562],
         [0, 0, 0, 1]] )}
];