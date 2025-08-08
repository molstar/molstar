/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { decodeColor } from '../../../extensions/mvs/helpers/utils';
import { MVSData_States } from '../../../extensions/mvs/mvs-data';
import { createMVSBuilder, Structure as MVSStructure, Representation, Root } from '../../../extensions/mvs/tree/mvs/mvs-builder';
import { MVSNodeParams } from '../../../extensions/mvs/tree/mvs/mvs-tree';
import { ColorT, ComponentExpressionT, isPrimitiveComponentExpressions, PrimitivePositionT } from '../../../extensions/mvs/tree/mvs/param-types';
import { Mat3, Mat4, Vec3 } from '../../../mol-math/linear-algebra';

const Domains = {
    ChainA: { auth_asym_id: 'A' },

    SH2: { auth_asym_id: 'A', beg_auth_seq_id: 146, end_auth_seq_id: 247 },
    SH3: { auth_asym_id: 'A', beg_auth_seq_id: 83, end_auth_seq_id: 145 },
    P_loop: { auth_asym_id: 'A', beg_auth_seq_id: 246, end_auth_seq_id: 255 },
    Activation_loop: { auth_asym_id: 'A', beg_auth_seq_id: 384, end_auth_seq_id: 402 },
};

const DomainColors = {
    SH2: '#8ED1A4' as ColorT,
    SH2_BCR: '#D03B4B' as ColorT,
    SH3: '#64B9AA' as ColorT,
    P_loop: 'pink' as ColorT,
    Activation_loop: 'red' as ColorT,
    DFG_motif: 'orange' as ColorT,
};

const Colors = {
    '1opl': '#4577B2' as ColorT,
    '2gqg': '#BC536D' as ColorT,
    '2g2i': '#BC536D' as ColorT,
    '1iep': '#B9E3A0' as ColorT,
    '3ik3': '#F3774B' as ColorT,
    '3oxz': '#7D7EA5' as ColorT,

    'active-site': '#F3794C' as ColorT,
    'binding-site': '#FEEB9F' as ColorT,
};

// Obtained using https://www.rcsb.org/alignment
// Aligned to 1iep
const Superpositions = {
    '1opl': [-0.6321036327, 0.3450463255, 0.6938213248, 0, -0.6288677634, -0.7515716885, -0.1991615756, 0, 0.4527364948, -0.5622126202, 0.6920597055, 0, 36.3924122492, 118.2516908402, -26.4992054179, 1] as unknown as Mat4,
    '3ik3': [-0.7767826245, -0.6295936551, 0.0148520572, 0, 0.6059737752, -0.7408035481, 0.2898376906, 0, -0.1714775143, 0.2341408391, 0.9569605684, 0, 21.0648276775, 53.0266628762, -0.3385906075, 1] as unknown as Mat4,
    '2gqg': [0.0648740828, -0.7163272638, 0.6947421137, 0, 0.0160329972, -0.6953706204, -0.7184724374, 0, 0.9977646498, 0.0577490387, -0.0336266582, 0, -31.0690973964, 146.0940883054, 39.7107422531, 1] as unknown as Mat4,
    '2g2i': [-0.5680242227, 0.6527660987, 0.5012433569, 0, -0.10067389, 0.5493518768, -0.8295042395, 0, -0.8168312251, -0.5216406194, -0.2463286704, 0, -8.1905690894, 75.7603329146, -6.1327389269, 1] as unknown as Mat4,
    '3oxz': [0.7989033646, 0.5984398921, -0.0601922711, 0, -0.1303123126, 0.269921501, 0.9540236289, 0, 0.5871729857, -0.754328893, 0.2936252816, 0, -8.0697093741, 58.1709160658, 19.0363028443, 1] as unknown as Mat4,
};

const Steps = [
    {
        header: 'Transition Test',
        key: 'transition_test',
        description: ``,
        linger_duration_ms: 2000,
        transition_duration_ms: 500,
        state: (): Root => {
            const builder = createMVSBuilder();

            const _1cbs = structure(builder, '1cbs');
            const [poly,] = polymer(_1cbs, { color: Colors['1opl'] });

            const surface = poly.representation({
                type: 'surface',
                surface_type: 'gaussian',
            })

            _1cbs.component({ selector: 'ligand' })
                .transform({ ref: 'xform', translation: [0, 20, 0] })
                .representation({ type: 'ball_and_stick' })
                .color({ color: 'red' });

            surface.clip({
                ref: 'clip',
                type: 'plane',
                point: [22.0, 15, 0],
                normal: [0, 0, 1],
            });

            builder.transition({
                ref: 'clip-transition',
                target_ref: 'clip',
                duration_ms: 2000,
                property: ['point', 2],
                target_value: 55,
            });

            builder.transition({
                // ref: 'clip-transition',
                target_ref: 'xform',
                duration_ms: 2000,
                property: ['translation', 1],
                target_value: 0,
            });

            return builder;
        },
        camera: undefined as MVSNodeParams<'camera'> | undefined
        // camera: {
        //     position: [103.72, 69.35, 20.52],
        //     target: [0.36, 55.32, 21.8],
        //     up: [-0.01, 0.01, -1],
        // } satisfies MVSNodeParams<'camera'>,
    },
    {
        header: 'Transition Test',
        key: 'transition_test2',
        description: ``,
        linger_duration_ms: 2000,
        transition_duration_ms: 0,
        state: (): Root => {
            const builder = createMVSBuilder();

            const _1cbs = structure(builder, '1cbs');
            const [poly,] = polymer(_1cbs, { color: Colors['1opl'] });

            poly.representation({
                type: 'surface',
                surface_type: 'gaussian',
            }).opacity({ ref: 'opacity', opacity: 1 });

            _1cbs.component({ selector: 'ligand' })
                .transform({ ref: 'xform', translation: [0, 0, 0] })
                .representation({ type: 'ball_and_stick' })
                .color({ color: 'red' });


            // surface.clip({
            //     ref: 'clip',
            //     type: 'plane',
            //     point: [22.0, 25, 55],
            //     normal: [0, 0, 1],
            // });

            builder.transition({
                target_ref: 'opacity',
                duration_ms: 1000,
                property: ['opacity'],
                target_value: 0.33,
            });

            return builder;
        },
        camera: undefined as MVSNodeParams<'camera'> | undefined
        // camera: {
        //     position: [103.72, 69.35, 20.52],
        //     target: [0.36, 55.32, 21.8],
        //     up: [-0.01, 0.01, -1],
        // } satisfies MVSNodeParams<'camera'>,
    }
];


type Interaction = [label: string, polymer: PrimitivePositionT, ligand: PrimitivePositionT, options?: { skipResidue?: boolean }]

function drawInteractions(structure: MVSStructure, interactions: Interaction[]) {
    const primitives = structure.primitives();

    const interactingResidues: ComponentExpressionT[] = [];
    const addedResidues = new Set<string>();

    for (const [tooltip, a, b, options] of interactions) {
        primitives.tube({ start: a, end: b, color: '#4289B5', tooltip, radius: 0.1, dash_length: 0.1 });

        if (options?.skipResidue) continue;

        const expressions = isPrimitiveComponentExpressions(a) ? a.expressions! : [a as ComponentExpressionT];
        for (const _e of expressions) {
            const e = { ..._e };
            delete e.auth_atom_id;
            delete e.label_atom_id;

            const key = JSON.stringify(e);
            if (addedResidues.has(key)) continue;
            interactingResidues.push(e);
            addedResidues.add(key);
        }
    }

    if (interactingResidues.length === 0) return;

    structure
        .component({ selector: interactingResidues })
        .representation({ type: 'ball_and_stick' })
        .color({
            custom: {
                molstar_color_theme_name: 'element-symbol',
                molstar_color_theme_params: { carbonColor: { name: 'element-symbol', params: {} } },
            }
        });
}

function transform(structure: MVSStructure, id: keyof typeof Superpositions) {
    const rotation = Mat3.fromMat4(Mat3.zero(), Superpositions[id]);
    const translation = Mat4.getTranslation(Vec3.zero(), Superpositions[id]) as any;
    return structure.transform({ rotation, translation });
}

function wwPDBLink(id: string) {
    return `https://doi.org/10.2210/pdb${id.toLowerCase()}/pdb`;
}

function structure(builder: Root, id: string): MVSStructure {
    let ret = builder
        .download({ url: pdbUrl(id) })
        .parse({ format: 'bcif' })
        .modelStructure();

    if (id in Superpositions) {
        ret = transform(ret, id as any);
    }

    return ret;
}

function domains(structure: MVSStructure, reprensentation: Representation, domains: [selector: ComponentExpressionT, color: ColorT, label?: string, options?: { label_size?: number }][], options?: { label_size?: number }) {
    const hasLabels = domains.some(d => !!d[2]);
    const primitives = hasLabels ? structure.primitives() : undefined;

    for (const [selector, color, label, opts] of domains) {
        reprensentation.color({ selector, color });
        if (label) primitives!.label({ position: selector, text: label, label_color: color, label_size: opts?.label_size ?? options?.label_size ?? 1.5 });
    }
}

function polymer(structure: MVSStructure, options: { color: ColorT }) {
    const component = structure.component({ selector: { label_asym_id: 'A' } });
    const reprensentation = component.representation({ type: 'cartoon' });
    reprensentation.color({ color: options.color });
    return [component, reprensentation] as const;
}

function ligand(structure: MVSStructure, options: {
    selector: ComponentExpressionT | ComponentExpressionT[],
    label?: string,
    surface?: boolean,
    carbon_color?: ColorT,
    uniform_color?: ColorT,
    label_color?: ColorT,
    label_size?: number,
    opacity?: number,
    component_ref?: string,
}) {
    const comp = structure.component({ ref: options?.component_ref, selector: options.selector });
    const coloring = options.uniform_color
        ? { color: options.uniform_color }
        : {
            custom: {
                molstar_color_theme_name: 'element-symbol',
                molstar_color_theme_params: { carbonColor: options?.carbon_color ? { name: 'uniform', params: { value: decodeColor(options?.carbon_color) } } : { name: 'element-symbol', params: {} } }
            }
        };

    if (options.surface) comp.representation({ type: 'surface' }).color(coloring).opacity({ opacity: 0.33 });
    const repr = comp.representation({ type: 'ball_and_stick' }).color(coloring);
    if (options.opacity) repr.opacity({ opacity: options.opacity });

    const label_color: ColorT = options?.label_color ?? options.uniform_color ?? options.carbon_color ?? '#5B53A4';
    if (options.label) {
        structure.primitives().label({
            position: Array.isArray(options.selector) ? { expressions: options.selector } : options.selector,
            text: options.label,
            label_color,
            label_size: options?.label_size ?? 1.5
        });
    }

    return comp;
}

function bindingSite(structure: MVSStructure, residues: [selector: ComponentExpressionT, label: string][], options: {
    color?: ColorT,
    label_size?: number,
}) {
    const color: ColorT = options.color ?? '#5B53A4';
    const coloring = {
        custom: {
            molstar_color_theme_name: 'element-symbol',
            molstar_color_theme_params: { carbonColor: { name: 'uniform', params: { value: decodeColor(color) } } }
        }
    };

    structure.component({ selector: residues.map(r => r[0]) }).representation({ type: 'ball_and_stick' }).color(coloring);

    const primitives = structure.primitives();
    for (const [selector, label] of residues) {
        primitives.label({
            position: selector,
            text: label,
            label_color: color,
            label_size: options?.label_size ?? 1.5
        });
    }
}

function pdbUrl(id: string) {
    return `https://www.ebi.ac.uk/pdbe/entry-files/download/${id.toLowerCase()}.bcif`;
}

export function buildStory(): MVSData_States {
    const snapshots = Steps.map((s, i) => {
        const builder = s.state();
        if (s.camera) builder.camera(s.camera);

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
            title: 'The Structural Story of BCR-ABL: A Kinase Out of Control',
            version: '1.0',
            timestamp: new Date().toISOString(),
        }
    };
}