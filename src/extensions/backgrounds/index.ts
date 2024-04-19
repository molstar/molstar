/**
 * Copyright (c) 2022-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginBehavior } from '../../mol-plugin/behavior/behavior';
import { PluginConfig } from '../../mol-plugin/config';
import { Color } from '../../mol-util/color/color';

// from https://visualsonline.cancer.gov/details.cfm?imageid=2304, public domain
import image_cells from './images/cells.jpg';

// created with http://alexcpeterson.com/spacescape/
import face_nebula_nx from './skyboxes/nebula/nebula_left2.jpg';
import face_nebula_ny from './skyboxes/nebula/nebula_bottom4.jpg';
import face_nebula_nz from './skyboxes/nebula/nebula_back6.jpg';
import face_nebula_px from './skyboxes/nebula/nebula_right1.jpg';
import face_nebula_py from './skyboxes/nebula/nebula_top3.jpg';
import face_nebula_pz from './skyboxes/nebula/nebula_front5.jpg';

export const Backgrounds = PluginBehavior.create<{ }>({
    name: 'extension-backgrounds',
    category: 'misc',
    display: {
        name: 'Backgrounds'
    },
    ctor: class extends PluginBehavior.Handler<{ }> {
        register(): void {
            this.ctx.config.set(PluginConfig.Background.Styles, [
                [{
                    variant: {
                        name: 'off',
                        params: {}
                    }
                }, 'Off'],
                [{
                    variant: {
                        name: 'radialGradient',
                        params: {
                            centerColor: Color(0xFFFFFF),
                            edgeColor: Color(0x808080),
                            ratio: 0.2,
                            coverage: 'viewport',
                        }
                    }
                }, 'Light Radial Gradient'],
                [{
                    variant: {
                        name: 'image',
                        params: {
                            source: {
                                name: 'url',
                                params: image_cells
                            },
                            lightness: 0,
                            saturation: 0,
                            opacity: 1,
                            blur: 0,
                            coverage: 'viewport',
                        }
                    }
                }, 'Normal Cells Image'],
                [{
                    variant: {
                        name: 'skybox',
                        params: {
                            faces: {
                                name: 'urls',
                                params: {
                                    nx: face_nebula_nx,
                                    ny: face_nebula_ny,
                                    nz: face_nebula_nz,
                                    px: face_nebula_px,
                                    py: face_nebula_py,
                                    pz: face_nebula_pz,
                                }
                            },
                            lightness: 0,
                            saturation: 0,
                            opacity: 1,
                            blur: 0.3,
                            rotation: { x: 0, y: 0, z: 0 },
                        }
                    }
                }, 'Purple Nebula Skybox'],
            ]);
        }

        update() {
            return false;
        }

        unregister() {
            this.ctx.config.set(PluginConfig.Background.Styles, []);
        }
    },
    params: () => ({ })
});
