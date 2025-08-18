/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MVSData_States } from '../../../extensions/mvs/mvs-data';
import { createMVSBuilder, Structure as MVSStructure, Root } from '../../../extensions/mvs/tree/mvs/mvs-builder';
import { MVSNodeParams } from '../../../extensions/mvs/tree/mvs/mvs-tree';
import { ColorT } from '../../../extensions/mvs/tree/mvs/param-types';
import { Mat4 } from '../../../mol-math/linear-algebra';

const Colors = {
    '1cbs': '#4577B2' as ColorT,

    'ligand-away': '#F3794C' as ColorT,
    'ligand-docked': '#B9E3A0' as ColorT,
};

const Steps = [
    {
        header: 'Animation Demo',
        key: 'intro',
        description: `### Molecular Animation
A story showcasing MolViewSpec animation capabilities.`,
        linger_duration_ms: 2000,
        transition_duration_ms: 500,
        state: (): Root => {
            const builder = createMVSBuilder();

            const _1cbs = structure(builder, '1cbs');
            polymer(_1cbs, { color: Colors['1cbs'] });

            const prims = _1cbs.primitives({
                ref: 'prims',
                label_opacity: 0,
            });
            prims.label({ text: 'Animation Demo', position: { label_asym_id: 'A' }, label_size: 10 });

            const anim = builder.animation({
                custom: {
                    molstar_trackball: {
                        name: 'rock',
                        params: { speed: 1.5 },
                    }
                }
            });
            anim.interpolate({
                kind: 'scalar',
                ref: 'prims-opacity',
                target_ref: 'prims',
                duration_ms: 1000,
                property: 'label_opacity',
                end: 1,
            });


            // Uncomment this to make 2nd frame render much faster
            // It will cause shader compilation to happen during the 1st snapshot

            // const surface = poly.representation({
            //     type: 'surface',
            //     surface_type: 'gaussian',
            // }).opacity({ opacity: 0 });

            // _1cbs.component({ selector: 'ligand' })
            //     .representation({ type: 'ball_and_stick' })
            //     .opacity({ opacity: 0 });

            // surface.clip({
            //     ref: 'clip',
            //     type: 'plane',
            //     point: [22.0, 15, 0],
            //     normal: [0, 0, 1],
            // });

            return builder;
        },
        camera: {
            position: [-11.49, -37.05, 15.78],
            target: [15.85, 17.26, 24.32],
            up: [-0.88, 0.4, 0.26],
        } satisfies MVSNodeParams<'camera'>,
    },
    {
        header: 'Ligand Docking',
        description: `Animate ligand moving to the binding site`,
        linger_duration_ms: 2500,
        transition_duration_ms: 500,
        state: (): Root => {
            const builder = createMVSBuilder();

            const _1cbs = structure(builder, '1cbs');
            const [poly,] = polymer(_1cbs, { color: Colors['1cbs'] });

            const surface = poly.representation({
                type: 'surface',
                surface_type: 'gaussian',
            });

            _1cbs.component({ selector: 'ligand' })
                .transform({
                    ref: 'xform',
                    translation: [5, 20, -20],
                    rotation: [1, 0, 0, 0, 1, 0, 0, 0, 1],
                    rotation_center: 'centroid',
                })
                .representation({ type: 'ball_and_stick' })
                .color({ ref: 'ligand-color', color: 'red' });

            surface.clip({
                ref: 'clip',
                type: 'plane',
                point: [22.0, 15, 0],
                normal: [0, 0, 1],
            });

            const anim = builder.animation();

            anim.interpolate({
                kind: 'scalar',
                ref: 'clip-transition',
                target_ref: 'clip',
                duration_ms: 2000,
                property: ['point', 2],
                end: 55,
                easing: 'sin-in',
            });

            anim.interpolate({
                kind: 'vec3',
                target_ref: 'xform',
                duration_ms: 2000,
                property: 'translation',
                end: [0, 0, 0],
                noise_magnitude: 1,
            });

            anim.interpolate({
                kind: 'rotation_matrix',
                target_ref: 'xform',
                duration_ms: 2000,
                property: 'rotation',
                noise_magnitude: 0.2,
            });

            anim.interpolate({
                kind: 'color',
                target_ref: 'ligand-color',
                duration_ms: 2000,
                property: 'color',
                end: Colors['ligand-docked'],
            });

            return builder;
        },
        camera: {
            position: [-30.63, 77.29, 2.28],
            target: [19.16, 26.15, 22.82],
            up: [0.69, 0.71, 0.09],
        } satisfies MVSNodeParams<'camera'>,
    },
    {
        header: 'Highlight & Opacity',
        description: `Animate emissive, opacity and transform properties`,
        linger_duration_ms: 2000,
        transition_duration_ms: 0,
        state: (): Root => {
            const builder = createMVSBuilder();

            const _1cbs = structure(builder, '1cbs');
            const [poly,] = polymer(_1cbs, { color: Colors['1cbs'] });

            poly.representation({
                type: 'surface',
                surface_type: 'gaussian'
            }).opacity({ ref: 'opacity', opacity: 1 }).color({ ref: 'surface-color', color: 'white' });

            _1cbs.component({ selector: 'ligand' })
                .transform({ ref: 'xform', translation: [0, 0, 0] })
                .representation({
                    ref: 'repr',
                    type: 'ball_and_stick',
                    custom: {
                        molstar_reprepresentation_params: {
                            emissive: 0,
                        }
                    }
                })
                .color({ color: Colors['ligand-docked'] });

            const primitives = builder.primitives({
                ref: 'primitives',
                instances: [
                    Mat4.identity()
                ],
                opacity: 0,
            });

            primitives.ellipsoid({
                center: [0, 0, 0],
                radius: [2, 3, 2.5],
                color: 'red'
            });

            const anim = builder.animation();

            anim.interpolate({
                kind: 'scalar',
                target_ref: 'repr',
                duration_ms: 1000,
                property: ['custom', 'molstar_reprepresentation_params', 'emissive'],
                end: 0.2,
            });

            anim.interpolate({
                kind: 'scalar',
                target_ref: 'opacity',
                duration_ms: 1000,
                frequency: 2,
                alternate_direction: true,
                property: 'opacity',
                end: 0,
            });

            anim.interpolate({
                kind: 'transform_matrix',
                target_ref: 'primitives',
                property: ['instances', 0],
                translation_start: [20.24, 29.64, 14.85],
                translation_end: [21.84, 21.71, 27.04],
                translation_frequency: 4,
                pivot: [0, 0, 0],
                rotation_noise_magnitude: 0.2,
                scale_end: [0.01, 0.01, 0.01],
                duration_ms: 1000,
            });

            anim.interpolate({
                kind: 'scalar',
                target_ref: 'primitives',
                duration_ms: 1000,
                property: 'opacity',
                end: 1,
            });


            anim.interpolate({
                kind: 'color',
                target_ref: 'surface-color',
                duration_ms: 2000,
                property: 'color',
                palette: {
                    kind: 'continuous',
                    colors: ['white', Colors['1cbs'], 'white'],
                }
            });

            return builder;
        },
        camera: {
            position: [6.92, 47.17, 10.68],
            target: [21.79, 22.2, 23.43],
            up: [0.8, 0.57, 0.2],
        } satisfies MVSNodeParams<'camera'>,
    }
];

function structure(builder: Root, id: string): MVSStructure {
    return builder
        .download({ url: pdbUrl(id) })
        .parse({ format: 'bcif' })
        .modelStructure();
}

function polymer(structure: MVSStructure, options: { color: ColorT }) {
    const component = structure.component({ selector: { label_asym_id: 'A' } });
    const reprensentation = component.representation({ type: 'cartoon' });
    reprensentation.color({ color: options.color });
    return [component, reprensentation] as const;
}

function pdbUrl(id: string) {
    return `https://www.ebi.ac.uk/pdbe/entry-files/download/${id.toLowerCase()}.bcif`;
}

export function buildStory(): MVSData_States {
    const snapshots = Steps.map((s, i) => {
        const builder = s.state();
        if (s.camera) builder.camera(s.camera);

        builder.canvas({
            custom: {
                molstar_postprocessing: {
                    enable_outline: true,
                    enable_ssao: true,
                }
            }
        });

        const description = i > 0 ? `${s.description}\n\n[Go to start](#intro)` : s.description;

        return builder.getSnapshot({
            title: s.header,
            key: s.key,
            description,
            description_format: 'markdown',
            linger_duration_ms: s.linger_duration_ms ?? 500,
            transition_duration_ms: s.transition_duration_ms ?? 1000,
        });
    });

    return {
        kind: 'multiple',
        snapshots,
        metadata: {
            title: 'Animation Showcase',
            version: '1.0',
            timestamp: new Date().toISOString(),
        }
    };
}