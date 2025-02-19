/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MVSData_States } from '../../../extensions/mvs/mvs-data';
import { createMVSBuilder, Root } from '../../../extensions/mvs/tree/mvs/mvs-builder';
import { Mat3, Mat4, Vec3 } from '../../../mol-math/linear-algebra';

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

            const _1opl = builder
                .download({ url: pdbUrl('1opl') })
                .parse({ format: 'bcif' })
                .modelStructure();

            const comp = _1opl.component({ selector: 'polymer' });
            comp.representation().color({ color: 'blue' });
            comp.label({ text: 'ABL Kinase' });

            return builder;
        },
    }, {
        header: 'The Birth of a Rogue Kinase',
        description: `
### The ABL Kinase: A Well-Regulated Enzyme

Normally, the ABL kinase ([PDB ID 1OPL](https://www.ebi.ac.uk/pdbe/entry/pdb/1opl/index)) is a well-regulated enzyme, kept in check by its SH3 and SH2 domains, which fold back onto the kinase domain like a safety lock.
`,
        state: () => {
            const builder = createMVSBuilder();
            const _1opl = builder
                .download({ url: pdbUrl('1opl') })
                .parse({ format: 'bcif' })
                .modelStructure();

            const SH3 = { auth_asym_id: 'A', beg_auth_seq_id: 83, end_auth_seq_id: 145 };
            const SH2 = { auth_asym_id: 'A', beg_auth_seq_id: 146, end_auth_seq_id: 247 };
            const NTerm = { auth_asym_id: 'A', beg_auth_seq_id: 259, end_auth_seq_id: 344 };
            const CTerm = { auth_asym_id: 'A', beg_auth_seq_id: 345, end_auth_seq_id: 519 };

            const polymer_1opl = _1opl.component({ selector: 'polymer' }).representation().color({ color: 'blue' });
            polymer_1opl.color({
                color: 'teal',
                selector: SH3,
            }).color({
                color: 'lightblue',
                selector: SH2,
            })
            // .color({
            //     color: 'gold',
            //     selector: NTerm,
            // }).color({
            //     color: 'pink',
            //     selector: CTerm,
            // });

            _1opl.component({ selector: SH3 }).label({ text: 'SH3 Domain' });
            _1opl.component({ selector: SH2 }).label({ text: 'SH2 Domain' });
            // _1opl.component({ selector: NTerm }).label({ text: 'N-terminal Lobe' });
            // _1opl.component({ selector: CTerm }).label({ text: 'C-terminal Lobe' });

            // _1opl.component({ selector: 'ligand' }).representation({ custom: { molstar_use_default_coloring: true } });

            _1opl.component({ selector: [SH2, SH3] }).focus();

            return builder;
        },
    }, {
        header: 'The Birth of a Rogue Kinase',
        description: `
### The Birth of a Rogue Kinase

But in BCR-ABL, this safety mechanism is gone. Comparing the full protein (in blue) to the kinase domain alone ([PDB ID 2GQG](https://www.ebi.ac.uk/pdbe/entry/pdb/2gqg/index), you can
see how the SH3 and SH2 domains (blue in normal ABL, red in BCR-ABL) are no longer positioned to restrain the kinase.

With this lock removed, BCR-ABL is stuck in an active conformation, like an accelerator pedal jammed to the floor. Without
its normal regulation, BCR-ABL will keep signaling, unchecked causing unregulated cell growth and cancer - chronic myeloid leukemia (CML).
`,
        state: () => {
            const builder = createMVSBuilder();

            const xform = [
                0.1786429439306385,
                0.36591659669680254,
                0.9133408117423878,
                0,
                -0.7115649550169003,
                0.689150100470433,
                -0.13692092093997132,
                0,
                -0.6795305062498747,
                -0.6254414985121173,
                0.3834853250021335,
                0,
                98.64312014906896,
                3.961551591952407,
                2.8217840460405625,
                1
            ];
            const rotation = Mat3.fromMat4(Mat3.zero(), xform as Mat4);
            const translation = Mat4.getTranslation(Vec3.zero(), xform as Mat4) as any;


            const _2gqg = builder
                .download({ url: pdbUrl('2gqg') })
                .parse({ format: 'bcif' })
                .modelStructure()
                .transform({ rotation, translation });

            _2gqg.component({ selector: 'polymer' }).representation().color({ color: 'red' });

            const _1opl = builder
                .download({ url: pdbUrl('1opl') })
                .parse({ format: 'bcif' })
                .modelStructure();

            const SH3 = { auth_asym_id: 'A', beg_auth_seq_id: 83, end_auth_seq_id: 145 };
            const SH2 = { auth_asym_id: 'A', beg_auth_seq_id: 146, end_auth_seq_id: 247 };

            const polymer_1opl = _1opl.component({ selector: 'polymer' }).representation().color({ color: 'blue' });
            polymer_1opl.color({
                color: 'teal',
                selector: SH3,
            }).color({
                color: 'lightblue',
                selector: SH2,
            });

            _1opl.component({ selector: SH3 }).label({ text: 'SH3 Domain' });
            _1opl.component({ selector: SH2 }).label({ text: 'SH2 Domain' });

            _1opl.component({ selector: { auth_asym_id: 'A' } }).focus();

            return builder;
        },
    }, {
        header: 'ATP Binding and Unstoppable Signaling',
        description: `
### ATP Binding and Unstoppable Signaling

To function, every kinase needs ATP, and BCR-ABL is no exception. Here, you can see non-hydrolysable ATP analogue (AMP-PNP) (green) nestled in the active site,
perfectly positioned for its phosphate group to be transferred to a substrate. Look closely at the active site residues—Lys271, Glu286, and Asp381 (in orange).
They form a crucial network that stabilizes ATP and catalyzes phosphorylation, allowing BCR-ABL to continuously activate downstream signaling pathways.

This interaction is key to orienting ATP for catalysis. In a normal kinase, ATP binding is a carefully controlled step.
But in BCR-ABL, there’s no regulation. ATP binds, reactions happen, and the leukemia-driving signals keep firing.
`,
        state: () => {
            const builder = createMVSBuilder();

            const xform = [
                0.1786429439306385,
                0.36591659669680254,
                0.9133408117423878,
                0,
                -0.7115649550169003,
                0.689150100470433,
                -0.13692092093997132,
                0,
                -0.6795305062498747,
                -0.6254414985121173,
                0.3834853250021335,
                0,
                98.64312014906896,
                3.961551591952407,
                2.8217840460405625,
                1
            ];
            const rotation = Mat3.fromMat4(Mat3.zero(), xform as Mat4);
            const translation = Mat4.getTranslation(Vec3.zero(), xform as Mat4) as any;

            const _2gqg = builder
                .download({ url: pdbUrl('2gqg') })
                .parse({ format: 'bcif' })
                .modelStructure()
                .transform({ rotation, translation });

            _2gqg.component({ selector: 'polymer' }).representation();

            const bindingSite = [
                { auth_asym_id: 'A', auth_seq_id: 271 },
                { auth_asym_id: 'A', auth_seq_id: 286 },
                { auth_asym_id: 'A', auth_seq_id: 381 },
            ];

            const ligand = { auth_asym_id: 'A', auth_seq_id: 501 };

            _2gqg.component({ selector: bindingSite }).representation({ type: 'ball_and_stick' }).color({ color: 'orange' });
            _2gqg.component({ selector: ligand }).representation({ type: 'ball_and_stick' }).color({ color: 'green' });

            const primitives = _2gqg.primitives();
            primitives.label({ position: bindingSite[0], text: 'Lys271', label_color: 'orange' });
            primitives.label({ position: bindingSite[1], text: 'Glu286', label_color: 'orange' });
            primitives.label({ position: bindingSite[2], text: 'Asp281', label_color: 'orange' });

            primitives.label({ position: ligand, text: 'AMP-PNP ATP analogue', label_color: 'green' });

            _2gqg.component({ selector: [...bindingSite, ligand] }).focus();

            return builder;
        },
    }, {
        header: 'Imatinib: The Drug That Changed Everything',
        description: `
### Imatinib: The Drug That Changed Everything

For years, chronic myeloid leukemia (CML) was a death sentence. Then came Imatinib (Gleevec), a molecule designed to fit into the ATP-binding pocket and lock BCR-ABL in
an inactive conformation. Viewing the Imatinib-bound structure (1IEP), and you'll see the difference—this time, the kinase is frozen. The drug (colour) forms key interactions
with Thr315, the gatekeeper residue; Asp381 and Phe382.

Notice how the P-loop (hot pink) (Viz 3.b), which normally cradles ATP, has shifted into a closed conformation and the activation loop is also flipped into its closed
conformation (green). Imatinib doesn’t just block ATP—it forces the kinase into a state where it can't function at all. Switching between the ATP-bound (active) and
Imatinib-bound (inactive) structures the change can be seen. The change is decisive: BCR-ABL is finally silenced.
`,
        state: () => {
            const builder = createMVSBuilder();

            const _1iep = builder
                .download({ url: pdbUrl('1iep') })
                .parse({ format: 'bcif' })
                .modelStructure();

            _1iep.component({ selector: 'polymer' }).representation();

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

            _1iep.component({ selector: gateway }).representation({ type: 'ball_and_stick' }).color({ color: 'orange' });
            _1iep.component({ selector: bindingSite }).representation({ type: 'ball_and_stick' }).color({ color: 'teal' });
            _1iep.component({ selector: ligand }).representation({ type: 'ball_and_stick' }).color({ color: 'green' });

            const primitives = _1iep.primitives();
            primitives.label({ position: gateway[0], text: 'Lys271', label_color: 'orange' });
            primitives.label({ position: gateway[1], text: 'Glu286', label_color: 'orange' });
            primitives.label({ position: bindingSite[0], text: 'Thr315', label_color: 'teal' });
            primitives.label({ position: bindingSite[1], text: 'Asp381', label_color: 'teal' });
            primitives.label({ position: bindingSite[2], text: 'Phe382', label_color: 'teal' });
            primitives.label({ position: ligand, text: 'Imatinib', label_color: 'green' });

            _1iep.component({ selector: [...bindingSite, ligand] }).focus();

            return builder;
        },
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
        }
    }
];

function pdbUrl(id: string) {
    return `https://www.ebi.ac.uk/pdbe/entry-files/download/${id.toLowerCase()}.bcif`;
}

export function buildStory(): MVSData_States {
    const snapshots = Steps.map(s => s.state().getSnapshot({
        title: s.header,
        key: s.key,
        description: s.description,
        description_format: 'markdown',
        linger_duration_ms: 5000,
    }));

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