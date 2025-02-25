/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MVSData_States } from '../../../extensions/mvs/mvs-data';
import { createMVSBuilder, Structure as MVSStructure, Root } from '../../../extensions/mvs/tree/mvs/mvs-builder';
import { MVSNodeParams } from '../../../extensions/mvs/tree/mvs/mvs-tree';
import { ColorT, ComponentExpressionT, isPrimitiveComponentExpressions, PrimitiveComponentExpressionT, PrimitivePositionT } from '../../../extensions/mvs/tree/mvs/param-types';
import { Mat3, Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { Color } from '../../../mol-util/color';


const Domains = {
    SH2_1opl: { auth_asym_id: 'A', beg_auth_seq_id: 146, end_auth_seq_id: 247 },
    SH3_1opl: { auth_asym_id: 'A', beg_auth_seq_id: 83, end_auth_seq_id: 145 },
    NTerm_1opl: { auth_asym_id: 'A', beg_auth_seq_id: 259, end_auth_seq_id: 344 },
    CTerm_1opl: { auth_asym_id: 'A', beg_auth_seq_id: 345, end_auth_seq_id: 519 },
};

const Imatinib_1iep = {
    N_lobe: { auth_asym_id: 'A', beg_auth_seq_id: 225, end_auth_seq_id: 321 },
    C_lobe: { auth_asym_id: 'A', beg_auth_seq_id: 322, end_auth_seq_id: 498 },
    P_loop: { auth_asym_id: 'A', beg_auth_seq_id: 246, end_auth_seq_id: 255 },
    Action_loop: { auth_asym_id: 'A', beg_auth_seq_id: 384, end_auth_seq_id: 402 },
    DFG_motif: { auth_asym_id: 'A', beg_auth_seq_id: 381, end_auth_seq_id: 383 },
};

const Colors = {
    SH2_1opl: '#A3DAA3' as ColorT,
    SH3_1opl: '#55ACAE' as ColorT,
    '1opl_A': '#4577B2' as ColorT,
    '1opl_B': '#438EB4' as ColorT,
    '2gqg_A': '#AB0F45' as ColorT,
    '2gqg_B': '#C02749' as ColorT,
    '2gqg-binding-site': '#F3794C' as ColorT,
    ligand: '#1A906D' as ColorT,
    background: '#A0A0A0' as ColorT,
    '1iep-binding-site': '#4595B4' as ColorT,
    '1iep-gateway': '#F3794C' as ColorT,
};

const TransformsTo1IEP = {
    '1opl': [-0.6321036327, 0.3450463255, 0.6938213248, 0, -0.6288677634, -0.7515716885, -0.1991615756, 0, 0.4527364948, -0.5622126202, 0.6920597055, 0, 36.3924122492, 118.2516908402, -26.4992054179, 1] as unknown as Mat4,
    '3ik3': [-0.7767826245, -0.6295936551, 0.0148520572, 0, 0.6059737752, -0.7408035481, 0.2898376906, 0, -0.1714775143, 0.2341408391, 0.9569605684, 0, 21.0648276775, 53.0266628762, -0.3385906075, 1] as unknown as Mat4,
    '2gqg': [0.0648740828, -0.7163272638, 0.6947421137, 0, 0.0160329972, -0.6953706204, -0.7184724374, 0, 0.9977646498, 0.0577490387, -0.0336266582, 0, -31.0690973964, 146.0940883054, 39.7107422531, 1] as unknown as Mat4,
    '3oxz': [0.7989033646, 0.5984398921, -0.0601922711, 0, -0.1303123126, 0.269921501, 0.9540236289, 0, 0.5871729857, -0.754328893, 0.2936252816, 0, -8.0697093741, 58.1709160658, 19.0363028443, 1] as unknown as Mat4,
};

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

function transform(structure: MVSStructure, id: keyof typeof TransformsTo1IEP) {
    const rotation = Mat3.fromMat4(Mat3.zero(), TransformsTo1IEP[id]);
    const translation = Mat4.getTranslation(Vec3.zero(), TransformsTo1IEP[id]) as any;
    return structure.transform({ rotation, translation });
}

function structure(builder: Root, id: string): MVSStructure {
    let ret = builder
        .download({ url: pdbUrl(id) })
        .parse({ format: 'bcif' })
        .modelStructure();

    if (id in TransformsTo1IEP) {
        ret = transform(ret, id as any);
    }

    return ret;
}

const Steps = [
    {
        header: 'A Kinase Out of Control',
        key: 'intro',
        description: `
### The Structural Story of BCR-ABL: A Kinase Out of Control

BCR-ABL is a classic case of how structural biology can drive drug discovery. By exploring these structures, you will understand:

- How a small genetic fusion creates a rogue kinase.
- How ATP binding fuels uncontrolled cancer growth.
- How Imatinib revolutionized treatment by locking the kinase in an inactive state.
- How a single mutation (T315I) enabled resistance and brought new challenges.
- How Ponatinib and future inhibitors are being designed to keep up in this ongoing battle.
`,
        state: (): Root => {
            const builder = createMVSBuilder();

            const _1opl = structure(builder, '1opl');
            const comp = _1opl.component({ selector: { auth_asym_id: 'A' } });
            comp
                .representation()
                .color({ selector: { auth_asym_id: 'A' }, color: Colors['1opl_A'] });
            comp.label({ text: 'ABL Kinase' });

            _1opl
                .component({ selector: { label_asym_id: 'D' } })
                .representation({ type: 'ball_and_stick' })
                .color({
                    custom: {
                        molstar_color_theme_name: 'element-symbol',
                        molstar_color_theme_params: { carbonColor: { name: 'uniform', params: { value: Color.fromHexStyle(Colors['1opl_A']) } } }
                    }
                });

            _1opl
                .component({ selector: { label_asym_id: 'C' } })
                .representation({ type: 'ball_and_stick' })
                .color({ color: Colors['1opl_A'] });

            return builder;
        },
        camera: {
            position: [167.49, 96.41, 33.61],
            target: [-2.17, 75.1, 34.08],
            up: [-0.01, 0.02, -1],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'The Birth of a Rogue Kinase [1/2]',
        description: `
### The ABL Kinase: A Well-Regulated Enzyme

Normally, the ABL kinase ([PDB ID 1OPL](https://www.ebi.ac.uk/pdbe/entry/pdb/1opl/index)) is a well-regulated enzyme, kept in check by its SH3 and SH2 domains, which fold back onto the kinase domain like a safety lock.
`,
        state: () => {
            const builder = createMVSBuilder();
            const _1opl = structure(builder, '1opl');

            _1opl
                .component({ selector: { auth_asym_id: 'A' } })
                .representation()
                .color({ selector: { auth_asym_id: 'A' }, color: Colors['1opl_A'] })
                // .color({ selector: { auth_asym_id: 'B' }, color: Colors['1opl_B'] })
                .color({ color: Colors.SH3_1opl, selector: Domains.SH3_1opl })
                .color({ color: Colors.SH2_1opl, selector: Domains.SH2_1opl });

            _1opl.component({ selector: Domains.SH3_1opl }).label({ text: 'SH3' });
            _1opl.component({ selector: Domains.SH2_1opl }).label({ text: 'SH2' });

            _1opl
                .component({ selector: { label_asym_id: 'D' } })
                .representation({ type: 'ball_and_stick' })
                .color({
                    custom: {
                        molstar_color_theme_name: 'element-symbol',
                        molstar_color_theme_params: { carbonColor: { name: 'uniform', params: { value: Color.fromHexStyle(Colors['1opl_A']) } } }
                    }
                });

            _1opl
                .component({ selector: { label_asym_id: 'C' } })
                .representation({ type: 'ball_and_stick' })
                .color({ color: Colors['1opl_A'] });

            return builder;
        },
        camera: {
            position: [83.27, -26.31, -3.54],
            target: [-0.93, 48.72, 12.58],
            up: [-0.33, -0.17, -0.93],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'The Birth of a Rogue Kinase [2/2]',
        description: `
### The Birth of a Rogue Kinase

But in BCR-ABL, this safety mechanism is gone. Comparing the full protein (in blue) to the kinase domain alone ([PDB ID 2GQG](https://www.ebi.ac.uk/pdbe/entry/pdb/2gqg/index)), you can
see how the SH3 and SH2 domains (blue in normal ABL, red in BCR-ABL) are no longer positioned to restrain the kinase.

With this lock removed, BCR-ABL is stuck in an active conformation, like an accelerator pedal jammed to the floor. Without
its normal regulation, BCR-ABL will keep signaling, unchecked causing unregulated cell growth and cancer - chronic myeloid leukemia (CML).
`,
        state: () => {
            const builder = createMVSBuilder();

            const _1opl = structure(builder, '1opl');

            _1opl
                .component({ selector: { auth_asym_id: 'A' } })
                .representation()
                .color({ selector: { auth_asym_id: 'A' }, color: Colors['1opl_A'] })
                // .color({ selector: { auth_asym_id: 'B' }, color: Colors['1opl_B'] })
                .color({ color: Colors.SH3_1opl, selector: Domains.SH3_1opl })
                .color({ color: Colors.SH2_1opl, selector: Domains.SH2_1opl });

            _1opl.component({ selector: Domains.SH3_1opl }).label({ text: 'SH3' });
            _1opl.component({ selector: Domains.SH2_1opl }).label({ text: 'SH2' });

            _1opl
                .component({ selector: { label_asym_id: 'D' } })
                .representation({ type: 'ball_and_stick' })
                .color({
                    custom: {
                        molstar_color_theme_name: 'element-symbol',
                        molstar_color_theme_params: { carbonColor: { name: 'uniform', params: { value: Color.fromHexStyle(Colors['1opl_A']) } } }
                    }
                });

            _1opl
                .component({ selector: { label_asym_id: 'C' } })
                .representation({ type: 'ball_and_stick' })
                .color({ color: Colors['1opl_A'] });

            const _2gqg = structure(builder, '2gqg');

            _2gqg
                .component({ selector: { auth_asym_id: 'A' } })
                .representation()
                .color({ selector: { auth_asym_id: 'A' }, color: Colors['2gqg_A'] });

            // .color({ selector: { auth_asym_id: 'B' }, color: Colors['2gqg_B'] });

            return builder;
        },
        camera: {
            position: [70.96, -42.72, 53.51],
            target: [15.94, 35.83, 21.87],
            up: [-0.06, -0.41, -0.91],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'ATP Binding and Unstoppable Signaling',
        description: `
### ATP Binding and Unstoppable Signaling

To function, every kinase needs ATP, and BCR-ABL is no exception. Here, you can see non-hydrolysable ATP analogue (AMP-PNP) (green) nestled in the active site,
perfectly positioned for its phosphate group to be transferred to a substrate. Look closely at the active site residues—Lys271, Glu286, and Asp381 (in orange).
They form a crucial network that stabilizes ATP and catalyzes phosphorylation, allowing BCR-ABL to continuously activate downstream signaling pathways.

This interaction is key to orienting ATP for catalysis. In a normal kinase, ATP binding is a carefully controlled step.
But in BCR-ABL, there's no regulation. ATP binds, reactions happen, and the leukemia-driving signals keep firing.
`,
        state: () => {
            const builder = createMVSBuilder();

            const _2gqg = structure(builder, '2gqg');

            _2gqg.component({ selector: { auth_asym_id: 'A' } }).representation().color({ color: Colors.background });

            const bindingSite = [
                { auth_asym_id: 'A', auth_seq_id: 271 },
                { auth_asym_id: 'A', auth_seq_id: 286 },
                { auth_asym_id: 'A', auth_seq_id: 381 },
            ];

            const ligand = { auth_asym_id: 'A', auth_seq_id: 501 };

            _2gqg.component({ selector: bindingSite }).representation({ type: 'ball_and_stick' }).color({ color: Colors['2gqg-binding-site'] });
            _2gqg.component({ selector: ligand }).representation({ type: 'ball_and_stick' }).color({ custom: { molstar_color_theme_name: 'element-symbol' } });

            const primitives = _2gqg.primitives();
            primitives.label({ position: bindingSite[0], text: 'Lys271', label_color: Colors['2gqg-binding-site'] });
            primitives.label({ position: bindingSite[1], text: 'Glu286', label_color: Colors['2gqg-binding-site'] });
            primitives.label({ position: bindingSite[2], text: 'Asp281', label_color: Colors['2gqg-binding-site'] });

            primitives.label({ position: ligand, text: 'AMP-PNP ATP analogue', label_color: Colors.ligand, label_size: 1.5 });

            return builder;
        },
        camera: {
            position: [46.94, 73.07, 15.53],
            target: [12.46, 60.95, 12.21],
            up: [0.18, -0.25, -0.95],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'Imatinib: The Drug That Changed Everything [1/2]',
        description: `
### Imatinib: The Drug That Changed Everything

For years, chronic myeloid leukemia (CML) was a death sentence. Then came Imatinib (Gleevec), a molecule designed to fit into the ATP-binding pocket
and lock BCR-ABL in an inactive conformation. Viewing the Imatinib-bound structure ([PDB ID 1IEP](https://www.ebi.ac.uk/pdbe/entry/pdb/1iep/index)), and you'll see the difference—this time, the kinase is frozen.
The drug (colour) forms key interactions with Thr315, the gatekeeper residue; Asp381 and Phe382.
`,
        state: () => {
            const builder = createMVSBuilder();

            const _1iep = structure(builder, '1iep');

            const _1iep_polymer = _1iep.component({ selector: { auth_asym_id: 'A' } });

            _1iep_polymer
                .representation()
                .color({ color: Colors.background });

            const gateway = [
                { auth_asym_id: 'A', auth_seq_id: 271 },
                { auth_asym_id: 'A', auth_seq_id: 286 },
            ];

            const bindingSite = [
                { auth_asym_id: 'A', auth_seq_id: 315 },
                { auth_asym_id: 'A', auth_seq_id: 381 },
                { auth_asym_id: 'A', auth_seq_id: 382 },
            ];

            const ligand = { auth_asym_id: 'A', auth_seq_id: 201 };

            _1iep.component({ selector: gateway }).representation({ type: 'ball_and_stick' }).color({ color: Colors['1iep-gateway'] });
            _1iep.component({ selector: bindingSite }).representation({ type: 'ball_and_stick' }).color({ color: Colors['1iep-binding-site'] });
            _1iep.component({ selector: ligand }).representation({ type: 'ball_and_stick' }).color({ custom: { molstar_color_theme_name: 'element-symbol' } });

            const primitives = _1iep.primitives();
            primitives.label({ position: gateway[0], text: 'Lys271', label_color: Colors['1iep-gateway'] });
            primitives.label({ position: gateway[1], text: 'Glu286', label_color: Colors['1iep-gateway'] });
            primitives.label({ position: bindingSite[0], text: 'Thr315', label_color: Colors['1iep-binding-site'] });
            primitives.label({ position: bindingSite[1], text: 'Asp381', label_color: Colors['1iep-binding-site'] });
            primitives.label({ position: bindingSite[2], text: 'Phe382', label_color: Colors['1iep-binding-site'] });
            primitives.label({ position: ligand, text: 'Imatinib', label_color: Colors.ligand, label_size: 1.5 });

            return builder;
        },
        camera: {
            position: [47.02, 69.08, 14.68],
            target: [12.92, 52.7, 17.97],
            up: [0.14, -0.46, -0.88],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'Imatinib: The Drug That Changed Everything [2/2]',
        description: `
### Imatinib: The Drug That Changed Everything

Notice how the P-loop (hot pink) (Viz 3.b), which normally cradles ATP, has shifted into a closed conformation and the activation loop is also flipped into its closed
conformation (green). Imatinib doesn't just block ATP—it forces the kinase into a state where it can't function at all. Switching between the ATP-bound (active) and
Imatinib-bound (inactive) structures the change can be seen. The change is decisive: BCR-ABL is finally silenced.
`,
        state: () => {
            const builder = createMVSBuilder();

            const _1iep = structure(builder, '1iep');

            const _1iep_polymer = _1iep.component({ selector: { auth_asym_id: 'A' } });

            _1iep_polymer
                .representation()
                .color({ color: Colors.background })
                .color({ color: 'lightgray', selector: Imatinib_1iep.N_lobe })
                .color({ color: 'darkgray', selector: Imatinib_1iep.C_lobe })
                .color({ color: 'pink', selector: Imatinib_1iep.P_loop })
                .color({ color: 'red', selector: Imatinib_1iep.Action_loop })
                .color({ color: 'orange', selector: Imatinib_1iep.DFG_motif });

            const _1iep_labels = _1iep.primitives();

            _1iep_labels.label({ position: Imatinib_1iep.P_loop, text: 'P-loop', label_size: 3, label_color: 'pink' });
            _1iep_labels.label({ position: Imatinib_1iep.N_lobe, text: 'N-lobe', label_size: 3, label_color: 'lightgray' });
            _1iep_labels.label({ position: Imatinib_1iep.C_lobe, text: 'C-lobe', label_size: 3, label_color: 'darkgray' });
            _1iep_labels.label({ position: Imatinib_1iep.Action_loop, text: 'Activation Loop', label_size: 3, label_color: 'red' });
            _1iep_labels.label({ position: Imatinib_1iep.DFG_motif, text: 'DFG Motif', label_size: 3, label_color: 'orange' });

            const gateway = [
                { auth_asym_id: 'A', auth_seq_id: 271 },
                { auth_asym_id: 'A', auth_seq_id: 286 },
            ];

            const bindingSite = [
                { auth_asym_id: 'A', auth_seq_id: 315 },
                { auth_asym_id: 'A', auth_seq_id: 381 },
                { auth_asym_id: 'A', auth_seq_id: 382 },
            ];

            const ligand = { auth_asym_id: 'A', auth_seq_id: 201 };

            _1iep.component({ selector: gateway }).representation({ type: 'ball_and_stick' }).color({ color: Colors['1iep-gateway'] });
            _1iep.component({ selector: bindingSite }).representation({ type: 'ball_and_stick' }).color({ color: Colors['1iep-binding-site'] });
            _1iep.component({ selector: ligand }).representation({ type: 'ball_and_stick' }).color({ custom: { molstar_color_theme_name: 'element-symbol' } });

            const primitives = _1iep.primitives();
            primitives.label({ position: gateway[0], text: 'Lys271', label_color: Colors['1iep-gateway'] });
            primitives.label({ position: gateway[1], text: 'Glu286', label_color: Colors['1iep-gateway'] });
            primitives.label({ position: bindingSite[0], text: 'Thr315', label_color: Colors['1iep-binding-site'] });
            primitives.label({ position: bindingSite[1], text: 'Asp381', label_color: Colors['1iep-binding-site'] });
            primitives.label({ position: bindingSite[2], text: 'Phe382', label_color: Colors['1iep-binding-site'] });
            primitives.label({ position: ligand, text: 'Imatinib', label_color: Colors.ligand, label_size: 1.5 });


            return builder;
        },
        camera: {
            position: [76.15, 47.15, 19.13],
            target: [12.92, 52.7, 17.97],
            up: [0.01, -0.08, -1],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'Resistance Strikes: The T315I Mutation [1/3]',
        description: `
### Resistance Strikes: The T315I Mutation

For a while, it seemed like leukemia had been beaten. But then, in some patients, the cancer returned. The culprit?
A single mutation: T315I, shown on [PDB ID 3IK3](https://www.ebi.ac.uk/pdbe/entry/pdb/3ik3/index) in orange. Here, you can see the mutation in red —what was once a threonine (Thr) is now an isoleucine (Ile).
It seems like a tiny change, but watch what happens when you compare this mutant structure to the Imatinib-bound one...
`,
        state: () => {
            const builder = createMVSBuilder();

            const mutation = { auth_asym_id: 'A', auth_seq_id: 315 };

            const _1iep = structure(builder, '1iep');

            _1iep
                .component({ selector: { auth_asym_id: 'A' } })
                .representation()
                .color({ color: '#BCE4A0' });

            _1iep.component({ selector: mutation }).representation({ type: 'ball_and_stick' }).color({ color: 'green' });

            const ligand = { auth_asym_id: 'A', auth_seq_id: 201 };
            _1iep.component({ selector: ligand }).representation({ type: 'ball_and_stick' }).color({ custom: { molstar_color_theme_name: 'element-symbol' } });

            const _3ik3 = structure(builder, '3ik3');

            _3ik3
                .component({ selector: { auth_asym_id: 'A' } })
                .representation()
                .color({ color: '#E75B49' });

            _3ik3.component({ selector: mutation }).representation({ type: 'ball_and_stick' }).color({ color: 'red' });
            const _3ik3_primitives = _3ik3.primitives();
            _3ik3_primitives.label({ position: mutation, text: 'T315I', label_color: 'red' });

            return builder;
        },
        camera: {
            position: [-14.69, 48.92, 15.26],
            target: [19.05, 55.11, 11.91],
            up: [-0.1, -0.02, -1],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'Resistance Strikes: The T315I Mutation [2/3]',
        description: `
### Resistance Strikes: The T315I Mutation

Thr315 was a crucial contact point for Imatinib. Now, with bulky isoleucine in its place, the drug's entry is blocked.
`,
        state: () => {
            const builder = createMVSBuilder();

            const mutation = { auth_asym_id: 'A', auth_seq_id: 315 };

            const _1iep = structure(builder, '1iep');

            _1iep
                .component({ selector: { auth_asym_id: 'A' } })
                .representation()
                .color({ color: '#BCE4A0' });

            _1iep.component({ selector: mutation }).representation({ type: 'ball_and_stick' }).color({ color: 'green' });

            const ligand = { auth_asym_id: 'A', auth_seq_id: 201 };
            _1iep.component({ selector: ligand }).representation({ type: 'spacefill' }).color({ custom: { molstar_color_theme_name: 'element-symbol' } });

            const _3ik3 = structure(builder, '3ik3');

            _3ik3
                .component({ selector: { auth_asym_id: 'A' } })
                .representation()
                .color({ color: '#E75B49' });

            _3ik3.component({ selector: mutation }).representation({ type: 'spacefill' }).color({ color: 'red' });

            return builder;
        },
        camera: {
            position: [-14.69, 48.92, 15.26],
            target: [19.05, 55.11, 11.91],
            up: [-0.1, -0.02, -1],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'Resistance Strikes: The T315I Mutation [3/3]',
        description: `
### Resistance Strikes: The T315I Mutation

The mutation is still allowing ATP to bind. The result? Resistance. BCR-ABL is active again, and the leukemia returns, this time untouchable by Imatinib.
`,
        state: () => {
            const builder = createMVSBuilder();

            const mutation = { auth_asym_id: 'A', auth_seq_id: 315 };

            const _2gqg = structure(builder, '2gqg');

            _2gqg.component({ selector: { auth_asym_id: 'A' } }).representation().color({ color: Colors.background });

            const ligand = { auth_asym_id: 'A', auth_seq_id: 501 };
            _2gqg.component({ selector: ligand }).representation({ type: 'spacefill' }).color({ custom: { molstar_color_theme_name: 'element-symbol' } });

            const primitives = _2gqg.primitives();
            primitives.label({ position: ligand, text: 'AMP-PNP ATP analogue', label_color: Colors.ligand, label_size: 1.5 });

            //

            const _3ik3 = structure(builder, '3ik3');

            _3ik3
                .component({ selector: { auth_asym_id: 'A' } })
                .representation()
                .color({ color: '#E75B49' });

            _3ik3.component({ selector: mutation }).representation({ type: 'spacefill' }).color({ color: 'red' });

            return builder;
        },
        camera: {
            position: [-14.69, 48.92, 15.26],
            target: [19.05, 55.11, 11.91],
            up: [-0.1, -0.02, -1],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'Fighting Back: Ponatinib and the Future of Kinase Inhibitors',
        description: `
### Fighting Back: Ponatinib and the Future of Kinase Inhibitors

The battle didn't end there. Scientists knew they needed a new inhibitor—one that could work even against T315I. Enter Ponatinib (shown in [PDB ID 3OXZ](https://www.ebi.ac.uk/pdbe/entry/pdb/3oxz/index)), a next-generation
drug designed to bypass resistance. Viewing the Ponatinib-bound structure, you'll see how it differs from Imatinib. Instead of being blocked by T315I,
Ponatinib has a flexible triple-bond linker, allowing it to slip into the binding site without clashing with the mutation.

Look closely at the interactions—Ponatinib forms new contacts that compensate for the loss of the Thr315 interaction. This structure tells a story of
rational drug design: scientists used everything they learned about BCR-ABL's structure to engineer a molecule that could fit where others failed.

But the story isn't over. New mutations continue to arise, and leukemia is still finding ways to outmaneuver our drugs. The future may lie in allosteric
inhibitors that bind outside the ATP pocket, or even in protein degradation strategies that eliminate BCR-ABL entirely. Whatever the next breakthrough is,
it will start here—with a deep understanding of structure and function, and the power of visualization to reveal the molecular battles happening
inside every cancer cell.
`,
        state: () => {
            const builder = createMVSBuilder();

            const mutation = { auth_asym_id: 'A', auth_seq_id: 315 };

            // const _1iep = structure(builder, '1iep');

            // _1iep
            //     .component({ selector: { auth_asym_id: 'A' } })
            //     .representation()
            //     .color({ color: '#BCE4A0' });

            // _1iep.component({ selector: mutation }).representation({ type: 'ball_and_stick' }).color({ color: 'green' });

            // const ligand = { auth_asym_id: 'A', auth_seq_id: 201 };
            // _1iep.component({ selector: ligand }).representation({ type: 'spacefill' }).color({ custom: { molstar_color_theme_name: 'element-symbol' } });

            const _3oxz = structure(builder, '3oxz');

            _3oxz
                .component({ selector: { auth_asym_id: 'A' } })
                .representation()
                .color({ color: Colors.background });

            _3oxz
                .component({ selector: mutation })
                .representation({ type: 'ball_and_stick' })
                .color({ color: 'red' });

            const _3ozx_ligand = { label_asym_id: 'B' };
            const _3oxz_lig = _3oxz
                .component({
                    selector: _3ozx_ligand,
                    // This is currently not supported in snapshots
                    // custom: { molstar_show_non_covalent_interactions: true, molstar_non_covalent_interactions_radius_ang: 5.0 },
                });

            const _3oxz_primitives = _3oxz.primitives();
            _3oxz_primitives.label({ position: mutation, text: 'T315I', label_color: 'red' });

            _3oxz_lig
                .representation({ type: 'ball_and_stick' })
                .color({ custom: { molstar_color_theme_name: 'element-symbol' } });

            drawInteractions(_3oxz, [
                ['H-bond', { auth_asym_id: 'A', auth_seq_id: 360, auth_atom_id: 'O' }, { label_asym_id: 'B', label_atom_id: 'N4' }],
                ['H-bond', { auth_asym_id: 'A', auth_seq_id: 361, auth_atom_id: 'O' }, { label_asym_id: 'B', label_atom_id: 'N4' }],
                ['H-bond', { auth_asym_id: 'A', auth_seq_id: 286, auth_atom_id: 'OE2' }, { label_asym_id: 'B', label_atom_id: 'N2' }],
                ['H-bond', { auth_asym_id: 'A', auth_seq_id: 381, auth_atom_id: 'N' }, { label_asym_id: 'B', label_atom_id: 'O1' }],
                ['H-bond', { auth_asym_id: 'A', auth_seq_id: 318, auth_atom_id: 'N' }, { label_asym_id: 'B', label_atom_id: 'N1' }],
                ['Pi-stacking',
                    { expressions: ['CD1', 'CD2', 'CE1', 'CE2', 'CG', 'CZ'].map(a => ({ auth_asym_id: 'A', auth_seq_id: 253, auth_atom_id: a })) },
                    { expressions: ['C81', 'C82', 'C83', 'C84', 'N81', 'N82'].map(a => ({ label_asym_id: 'B', auth_atom_id: a })) },
                ],
            ]);

            return builder;
        },
        camera: {
            position: [-9.53, 25.51, 22.64],
            target: [16.58, 52.22, 17.48],
            up: [0.35, -0.5, -0.79],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'The End',
        key: 'end',
        description: `
### The End

That's all folks! We hope you enjoyed this interactive journey through the structural biology of BCR-ABL.

The next time you look at a macromolecular structure, remember: each atom tells a story, and each discovery shapes the future of medicine.

Read more [here](https://pmc.ncbi.nlm.nih.gov/articles/PMC3513788/).
`,
        state: (): Root => {
            return Steps[0].state();
        },
        camera: {
            position: [167.49, 96.41, 33.61],
            target: [-2.17, 75.1, 34.08],
            up: [-0.01, 0.02, -1],
        } satisfies MVSNodeParams<'camera'>,
    }
];

function pdbUrl(id: string) {
    return `https://www.ebi.ac.uk/pdbe/entry-files/download/${id.toLowerCase()}.bcif`;
}

export function buildStory(): MVSData_States {
    const snapshots = Steps.map(s => {
        const builder = s.state();
        if (s.camera) builder.camera(s.camera);

        return builder.getSnapshot({
            title: s.header,
            key: s.key,
            description: s.description,
            description_format: 'markdown',
            linger_duration_ms: 5000,
            transition_duration_ms: 1500,
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