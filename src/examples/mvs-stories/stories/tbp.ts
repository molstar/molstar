/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { MVSData_States } from '../../../extensions/mvs/mvs-data';
import { createMVSBuilder, Structure as MVSStructure, Root } from '../../../extensions/mvs/tree/mvs/mvs-builder';
import { MVSNodeParams } from '../../../extensions/mvs/tree/mvs/mvs-tree';
import {
    ColorT,
    ComponentExpressionT,
    isPrimitiveComponentExpressions,
    PrimitivePositionT
} from '../../../extensions/mvs/tree/mvs/param-types';
import { Mat3, Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { decodeColor } from '../../../extensions/mvs/helpers/utils';

const Colors = {
    '1vok': '#4577B2' as ColorT,
    '1cdw': '#BC536D' as ColorT,
    '1cdw-2': '#c5a3af' as ColorT,
    '1vtl': '#B9E3A0' as ColorT,
    '7enc': '#0072B2' as ColorT,
    '7enc-2': '#D55E00' as ColorT,
    '7enc-3': '#009E73' as ColorT,
    '7enc-4': '#56B4E9' as ColorT,
};

// Obtained using https://www.rcsb.org/alignment
const Superpositions = {
    '1cdw': [-0.4665815186, 0.6063873444, -0.6438913535, 0, -0.581544075, -0.7588303199, -0.2932286385, 0, -0.6664144171, 0.2376361381, 0.7066971703, 0, 135.0863694935, 105.5007997009, 153.6890178993, 1] as unknown as Mat4,
    '1vtl': [-0.4769460004, 0.7214347188, -0.5020502557, 0, -0.297882932, 0.4047204695, 0.8645617968, 0, 0.8269149119, 0.5619014932, 0.0218732801, 0, 65.6043682658, -3.7328402905, -16.8650755387, 1] as unknown as Mat4,
    '7enc': [0.8975055044, -0.4316347566, -0.0904174009, 0, 0.247274877, 0.3227849997, 0.9136000105, 0, -0.3651561375, -0.8423189899, 0.3964337454, 0, -189.7572972798, 304.0841220076, -411.5005782853, 1] as unknown as Mat4,
};

const Steps = [
    {
        header: 'TATA-Binding Protein',
        key: 'intro',
        description: `
### TATA-Binding Protein Tells RNA Polymerase Where To Get Started on a Gene

Specialized DNA sequences next to genes, called promoters, define the proper start site and direction for transcription.
Promoters vary in sequence and location from organism to organism.
In eukaryotic cells, a complex promoter system that uses dozens of different proteins ensures that the proper RNA polymerase is targeted to each gene.
The TATA-binding protein (TBP) is the central element of this system, participating in transcription by all three RNA polymerases.

The crystal structure of TBP from *Arabidopsis thaliana* in its apo form is shown ([PDB ID 1VOK](${wwPDBLink('1vok')})).
The structure reveals a highly symmetric DNA-binding fold composed of two topologically identical domains derived from the two direct repeats in the phylogenetically conserved sequence.
The domains are related by an approximate 2-fold axis, and consist of a five-stranded, antiparallel beta sheet and two alpha helices.
It is thought that an ancient gene duplication created this protein by combining two copies of the same gene. 

The intramolecular symmetry generates a saddle-shaped structure dominated by a curved antiparallel beta sheet, which forms the saddle's concave face and interacts with the DNA.
The convex face of the saddle would interact with other proteins during transcription initiation.
`,
        state: (): Root => {
            const builder = createMVSBuilder();

            const _1vok = structure(builder, '1vok');
            const _1vok_comp = _1vok.component({ selector: { label_asym_id: 'A' } });
            _1vok_comp.representation({ type: 'cartoon' })
                .color({
                    custom: {
                        molstar_color_theme_name: 'sequence-id',
                        molstar_color_theme_params: { carbonColor: { name: 'sequence-id', params: {} } },
                    }
                });
            _1vok.primitives().label({
                position: { label_asym_id: 'A', label_seq_id: 88, label_atom_id: 'OD2' },
                text: 'TATA-Binding Protein',
                label_size: 5,
                label_color: Colors['1vok']
            });

            return builder;
        },
        camera: {
            position: [155.18, 118.49, 49.18],
            target: [74.81, 66.8, 30.7],
            up: [-0.55, 0.83, 0.1],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'TBP: Highly Conserved in Eukaryotes [1/2]',
        key: 'highly-conserved-1',
        description: `
### TATA-Binding Protein Is Highly Conserved in Eukaryotes

TBP has a phylogenetically conserved, 180 amino-acid carboxy-terminal domain (ranging from 38 to 93% identity among eukaryotes and archaebacteria), containing two structural repeats flanking a highly basic segment known as the basic repeat.
The C-terminal or core portion of the protein binds the TATA consensus sequence with high affinity, interacting with the minor groove and promoting DNA bending.
Structural superposition of TBP bound to DNA in human ([PDB ID 1CDW](${wwPDBLink('1cdw')})) and *A. thaliana* ([PDB ID 1VTL](${wwPDBLink('1vtl')})) shows that their sequences are 83% identical and the RMSD between the structures is 0.43 Å.
`,
        state: () => {
            const builder = createMVSBuilder();

            const _1cdw = structure(builder, '1cdw');
            select(_1cdw, { color: Colors['1cdw'], selector: 'protein' });
            select(_1cdw, { color: Colors['1cdw'], selector: 'nucleic', opacity: 0.5 })[0].label({ text: 'DNA' });
            // uses a range to 'adjust' font size of label
            label(_1cdw, { selector: { label_asym_id: 'C', beg_label_seq_id: 160, end_label_seq_id: 177 }, text: 'TBP (H. sapiens)' });

            const _1vtl = structure(builder, '1vtl');
            select(_1vtl, { color: Colors['1vtl'], selector: { label_asym_id: 'E' } });
            select(_1vtl, { color: Colors['1vtl'], selector: { label_asym_id: 'A' }, opacity: 0.5 });
            select(_1vtl, { color: Colors['1vtl'], selector: { label_asym_id: 'B' }, opacity: 0.5 });
            // uses a range to 'adjust' font size of label
            label(_1vtl, { selector: { label_asym_id: 'E', beg_label_seq_id: 75, end_label_seq_id: 92 }, text: 'TBP (A. thaliana)' });

            return builder;
        },
        camera: {
            position: [122.15, 104.37, 72.23],
            target: [77.48, 59.61, 30.36],
            up: [-0.73, 0.68, 0.06],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'TBP: Minor Groove [2/2]',
        key: 'highly-conserved-2',
        description: `
### TATA-Binding Protein Is Highly Conserved in Eukaryotes

Eukaryotic protein-coding genes transcribed by RNA polymerase II (pol II)  have a characteristic sequence of nucleotides, termed the TATA box, in front of the start site of transcription.
The typical sequence is something like T-A-T-A-a/t-A-a/t, where a/t refers to positions that can be either A or T.
TBP recognizes this TATA sequence and binds to it, creating a landmark that directs pol II to the transcription start site. 
The structures of *A. thaliana* and human core TBP-TATA element co-crystal structures demonstrate a common induced-fit mechanism of protein-DNA recognition involving subtle conformation changes in the protein and an unprecedented DNA distortion.
`,
        state: () => {
            const builder = createMVSBuilder();

            const _1cdw = structure(builder, '1cdw');
            select(_1cdw, { color: Colors['1cdw'], selector: 'protein' });
            select(_1cdw, { color: Colors['1cdw'], selector: 'nucleic', opacity: 0.5 });
            // select(_1cdw, { color: Colors['1cdw'], selector: { label_asym_id: 'A', beg_label_seq_id: 6, end_label_seq_id: 11 } })[0].label({ text: 'TATA box' });
            label(_1cdw, { selector: { label_asym_id: 'A', label_seq_id: 5 }, text: 'T' });
            label(_1cdw, { selector: { label_asym_id: 'A', label_seq_id: 6 }, text: 'A' });
            label(_1cdw, { selector: { label_asym_id: 'A', label_seq_id: 7 }, text: 'T' });
            label(_1cdw, { selector: { label_asym_id: 'A', label_seq_id: 8 }, text: 'A' });
            label(_1cdw, { selector: { label_asym_id: 'A', label_seq_id: 9 }, text: 'A' });
            label(_1cdw, { selector: { label_asym_id: 'A', label_seq_id: 10 }, text: 'A' });

            const _1vtl = structure(builder, '1vtl');
            select(_1vtl, { color: Colors['1vtl'], selector: { label_asym_id: 'E' } });
            select(_1vtl, { color: Colors['1vtl'], selector: { label_asym_id: 'A' }, opacity: 0.5 });
            select(_1vtl, { color: Colors['1vtl'], selector: { label_asym_id: 'B' }, opacity: 0.5 });

            return builder;
        },
        camera: {
            position: [108.48, 92.12, 4.1],
            target: [80.16, 56.15, 28.96],
            up: [-0.71, 0.68, 0.17],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'TBP: Binding to TATA Box [1/5]',
        key: 'tata-box-overview',
        description: `
### TATA-Binding Protein Binds to the Minor Groove of the TATA Box

When the first structures of TBP-DNA complexes were determined, researchers discovered that TBP is not gentle when it binds to DNA.
Instead, it grabs the TATA sequence, bends and unwinds it to open up the minor groove, and kinks it sharply in two places (*e.g.*, [PDB ID 1VTL](${wwPDBLink('1vtl')})) and ([PDB ID 1CDW](${wwPDBLink('1cdw')})). 

Interactions with the minor groove can be divided into different classes as seen in PDB ([PDB ID 1CDW](${wwPDBLink('1cdw')})).
The combination of interactions allows TATA-binding protein to recognize the proper DNA sequence.
`,
        state: () => {
            const builder = createMVSBuilder();

            const _1cdw = structure(builder, '1cdw');
            select(_1cdw, { color: Colors['1cdw'], selector: { label_asym_id: 'C' } });
            select(_1cdw, { color: Colors['1cdw-2'], selector: { label_asym_id: 'A' } });
            select(_1cdw, { color: Colors['1cdw-2'], selector: { label_asym_id: 'B' } });

            return builder;
        },
        camera: {
            position: [122.15, 104.37, 72.23],
            target: [77.48, 59.61, 30.36],
            up: [-0.73, 0.68, 0.06],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'TBP: Arginine [2/5]',
        key: 'tata-box-1',
        description: `
### Arginine Residues

A string of arginine amino acids interact with the phosphate groups of the DNA and glues the protein to the DNA: Arg192, Arg199, Arg204, and Arg290.
`,
        state: () => {
            const builder = createMVSBuilder();

            const _1cdw = structure(builder, '1cdw');
            select(_1cdw, { color: Colors['1cdw'], selector: { label_asym_id: 'C' } });
            select(_1cdw, { color: Colors['1cdw-2'], selector: { label_asym_id: 'A' } });
            select(_1cdw, { color: Colors['1cdw-2'], selector: { label_asym_id: 'B' } });

            label(_1cdw, { selector: { label_asym_id: 'C', label_seq_id: 38 }, text: 'Arg192' });
            drawInteractions(_1cdw, [
                ['H-bond', { label_asym_id: 'C', label_seq_id: 38, label_atom_id: 'NH2' }, { label_asym_id: 'B', label_seq_id: 7, label_atom_id: 'OP1' }],
                ['H-bond', { label_asym_id: 'C', label_seq_id: 38, label_atom_id: 'NH2' }, { label_asym_id: 'B', label_seq_id: 6, label_atom_id: `O3'` }],
                ['H-bond', { label_asym_id: 'C', label_seq_id: 38, label_atom_id: 'NH1' }, { label_asym_id: 'B', label_seq_id: 7, label_atom_id: 'OP1' }],
            ]);

            label(_1cdw, { selector: { label_asym_id: 'C', label_seq_id: 45 }, text: 'Arg199' });
            drawInteractions(_1cdw, [
                ['H-bond', { label_asym_id: 'C', label_seq_id: 45, label_atom_id: 'NE' }, { label_asym_id: 'B', label_seq_id: 8, label_atom_id: 'OP1' }],
                ['H-bond', { label_asym_id: 'C', label_seq_id: 45, label_atom_id: 'NH2' }, { label_asym_id: 'B', label_seq_id: 8, label_atom_id: 'OP1' }],
            ]);

            label(_1cdw, { selector: { label_asym_id: 'C', label_seq_id: 136 }, text: 'Arg290' });
            drawInteractions(_1cdw, [
                ['H-bond', { label_asym_id: 'C', label_seq_id: 136, label_atom_id: 'NH1' }, { label_asym_id: 'A', label_seq_id: 8, label_atom_id: 'OP1' }],
            ]);

            label(_1cdw, { selector: { label_asym_id: 'C', label_seq_id: 50 }, text: 'Arg204' });
            drawInteractions(_1cdw, [
                ['H-bond', { label_asym_id: 'B', label_seq_id: 9, label_atom_id: 'OP1' }, { label_asym_id: 'C', label_seq_id: 50, label_atom_id: 'NH1' }],
                ['H-bond', { label_asym_id: 'B', label_seq_id: 9, label_atom_id: 'OP1' }, { label_asym_id: 'C', label_seq_id: 50, label_atom_id: 'NH2' }],
            ]);

            return builder;
        },
        camera: {
            position: [113.87, 71.89, 26.29],
            target: [77.29, 61.61, 18.73],
            up: [-0.28, 0.96, 0.04],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'TBP: Phenylalanine [3/5]',
        key: 'tata-box-2',
        description: `
### Phenylalanine-Induced Kinks

At either end of the TATA element there are two pairs of phenylalanine side chains partially inserted between adjacent base pairs, producing the two dramatic kinks (Phe193, Phe210, Phe284, and Phe301). Between the two kinks the DNA is partially unwound.
`,
        state: () => {
            const builder = createMVSBuilder();

            const _1cdw = structure(builder, '1cdw');
            select(_1cdw, { color: Colors['1cdw'], selector: { label_asym_id: 'C' } });
            select(_1cdw, { color: Colors['1cdw-2'], selector: { label_asym_id: 'A' } });
            select(_1cdw, { color: Colors['1cdw-2'], selector: { label_asym_id: 'B' } });

            bindingSite(_1cdw, [
                [{ label_asym_id: 'C', label_seq_id: 39 }, 'Phe193'],
                [{ label_asym_id: 'C', label_seq_id: 56 }, 'Phe210'],
                [{ label_asym_id: 'C', label_seq_id: 130 }, 'Phe284'],
                [{ label_asym_id: 'C', label_seq_id: 147 }, 'Phe301'],
            ], { color: Colors['1cdw'] });

            return builder;
        },
        camera: {
            position: [111.01, 69.92, -3.14],
            target: [85.52, 57.53, 19.15],
            up: [-0.51, 0.85, -0.11],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'TBP: H-Bonds in Minor Groove [4/5]',
        key: 'tata-box-3',
        description: `
### Hydrogen Bonds in the Minor Groove

Polar side chains make minor groove hydrogen bonds with acceptors of base pairs centred about the approximate 2-fold symmetry axis (Asn163, Asn253, and Thr309).
`,
        state: () => {
            const builder = createMVSBuilder();

            const _1cdw = structure(builder, '1cdw');
            select(_1cdw, { color: Colors['1cdw'], selector: { label_asym_id: 'C' } });
            select(_1cdw, { color: Colors['1cdw-2'], selector: { label_asym_id: 'A' } });
            select(_1cdw, { color: Colors['1cdw-2'], selector: { label_asym_id: 'B' } });

            // Asn: 9, 99, 155
            label(_1cdw, { selector: { label_asym_id: 'C', label_seq_id: 9 }, text: 'Asn163' });
            label(_1cdw, { selector: { label_asym_id: 'C', label_seq_id: 155 }, text: 'Thr309' });
            drawInteractions(_1cdw, [
                ['H-bond', { label_asym_id: 'C', label_seq_id: 9, label_atom_id: 'ND2' }, { label_asym_id: 'B', label_seq_id: 8, label_atom_id: 'O2' }],
                ['H-bond', { label_asym_id: 'B', label_seq_id: 9, label_atom_id: 'O2' }, { label_asym_id: 'C', label_seq_id: 9, label_atom_id: 'ND2' }],
                ['H-bond', { label_asym_id: 'C', label_seq_id: 155, label_atom_id: 'OG1' }, { label_asym_id: 'A', label_seq_id: 8, label_atom_id: 'N3' }],
            ]);
            bindingSite(_1cdw, [
                [{ label_asym_id: 'C', label_seq_id: 99 }, 'Asn253'],
            ], { color: 'gray', label_color: Colors['1cdw'] });

            return builder;
        },
        camera: {
            position: [69.71, 56.65, 17.62],
            target: [75.68, 63.48, 32.29],
            up: [-0.9, 0.39, 0.18],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'TBP: Non-Polar Interactions [5/5]',
        key: 'tata-box-4',
        description: `
### Hydrophobic and van der Waals Interactions

Several residues projecting from the concave surface of TBP make hydrophobic or van der Waals side-chain/base contacts (<4 Å between non-hydrogen atoms) with the TATA element.
`,
        state: () => {
            const builder = createMVSBuilder();

            const _1cdw = structure(builder, '1cdw');
            select(_1cdw, { color: Colors['1cdw'], selector: { label_asym_id: 'C' } });
            select(_1cdw, { color: Colors['1cdw-2'], selector: { label_asym_id: 'A' } });
            select(_1cdw, { color: Colors['1cdw-2'], selector: { label_asym_id: 'B' } });

            surface(_1cdw, {
                selector: { label_asym_id: 'A' },
                carbon_color: Colors['1cdw-2'],
            });
            surface(_1cdw, {
                selector: { label_asym_id: 'B' },
                carbon_color: Colors['1cdw-2'],
            });

            return builder;
        },
        camera: {
            position: [122.15, 104.37, 72.23],
            target: [77.48, 59.61, 30.36],
            up: [-0.73, 0.68, 0.06],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'TBP and Transcription Pre-Initiation Complex',
        key: 'pre-init-complex',
        description: `
### TATA-Binding Protein and the Transcription Pre-Initiation Complex

The structure shown here, from [PDB ID 7ENC](${wwPDBLink('7enc')}), is that of a pol II “pre-initiation complex” (PIC) poised to start transcription.
This complex includes eukaryotic Mediator and pol II, along with a collection of general transcription factors that perform the necessary tasks of recognizing the transcription start site of genes, separating the two strands of the DNA double helix, and facilitating transcription initiation by RNA polymerase II.

TBP starts the process of assembling the PIC and works with TBP-associated factors (TAFs) as part of a large multi-component transcription factor, TFIID, that belongs to the collection of general transcription factors that cradle pol II with in the PIC.
After TBP binds to the promoter, it recruits additional transcription factors.
TFIIA and TFIIB interact with surrounding regions of the DNA and, along with TFIIF, assist with positioning of RNA polymerase II at the transcription start site.
TFIIE, TFIIH, and other TFIID subunits bring additional functionality to the complex.
In particular, part of TFIIH is a translocase that separates the two strands of DNA in preparation for transcription, and the CAK module of TFIIH adds phosphate groups to the long tail of RNA polymerase II, sending the signal that it is time to get started with mRNA synthesis. 

TBP is a critical player that gets the transcription process started and is central to the pre-initiation complex.
Lying at the very heart of this complex macromolecular machine, it is highly conserved among eukaryotes and *archaea*. 
`,
        state: () => {
            const builder = createMVSBuilder();

            const _7enc = structure(builder, '7enc');
            select(_7enc, { color: Colors['7enc'], selector: { label_asym_id: 'CB' } })[0].label({ text: 'TBP' });
            select(_7enc, { color: Colors['7enc-2'], selector: { label_asym_id: 'GB' } });
            select(_7enc, { color: Colors['7enc-2'], selector: { label_asym_id: 'HB' } });
            select(_7enc, { color: Colors['7enc-3'], opacity: 0.5 });

            _7enc.primitives()
                .label({ position: { label_entity_id: '57' }, text: 'pol II', label_size: 20, label_color: Colors['7enc'] })
                .label({ position: { label_entity_id: '53' }, text: 'DNA', label_size: 20, label_color: Colors['7enc-2'] })
                .label({ position: { label_entity_id: '21' }, text: 'mediator', label_size: 20, label_color: Colors['7enc-3'] })
                .label({ position: { label_entity_id: '42' }, text: 'TBP-associated factors (TAFs)', label_size: 20, label_color: Colors['7enc-4'] })
                .label({ position: [100, 150, 0], text: 'PIC', label_size: 30, label_color: Colors['1vok'] })
            ;

            return builder;
        },
        camera: {
            position: [122.15, 104.37, 72.23],
            target: [77.48, 59.61, 30.36],
            up: [-0.73, 0.68, 0.06],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'The End',
        key: 'end',
        description: `
### The End

That's all folks! We hope you enjoyed this interactive journey through the structural biology of the TATA-binding protein.

The next time you look at a macromolecular structure, remember: each atom tells a story, and each discovery shapes the future of medicine.

Read more in the relevant publications on PDB IDs [1VOK](https://doi.org/10.1038/nsb0994-621), [1CDW](https://doi.org/10.1073/pnas.93.10.4862), [1VTL](https://doi.org/10.1038/365520a0), and [7ENC](https://doi.org/10.1126/science.abg0635) as well as the Molecule of the Month articles on the [TATA-binding protein](https://pdb101.rcsb.org/motm/67) and [Mediator](https://pdb101.rcsb.org/motm/289).
`,
        state: (): Root => {
            return Steps[Steps.length - 2].state();
        },
        camera: {
            position: [461.95, 218.66, 140.67],
            target: [87.56, -22.14, 54.59],
            up: [-0.55, 0.83, 0.1],
        } satisfies MVSNodeParams<'camera'>,
    }
];

type Interaction = [label: string, res1: PrimitivePositionT, res2: PrimitivePositionT, options?: { skipResidue?: boolean }]

function drawInteractions(structure: MVSStructure, interactions: Interaction[]) {
    const primitives = structure.primitives();

    const interactingResidues: ComponentExpressionT[] = [];
    const addedResidues = new Set<string>();

    function drawResidue(a: PrimitivePositionT) {
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

    for (const [tooltip, a, b, options] of interactions) {
        primitives.tube({ start: a, end: b, color: '#4289B5', tooltip, radius: 0.1, dash_length: 0.1 });

        if (options?.skipResidue) continue;

        drawResidue(a);
        drawResidue(b);
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

function select(structure: MVSStructure, { color, opacity = 1.0, selector = 'polymer' }: { color: ColorT, opacity?: number, selector?: Partial<MVSNodeParams<'component'>>['selector'] }) {
    const component = structure.component({ selector });
    const representation = component.representation({ type: 'cartoon' });
    representation.color({ color }).opacity({ opacity });
    return [component, representation] as const;
}

function label(structure: MVSStructure, { selector, text }: { selector: Partial<MVSNodeParams<'component'>>['selector'], text: string }) {
    structure.component({ selector })
        .label({ text });
}

function surface(structure: MVSStructure, options: {
    selector: ComponentExpressionT | ComponentExpressionT[],
    carbon_color?: ColorT,
}) {
    const comp = structure.component({ selector: options.selector });
    const coloring = {
        custom: {
            molstar_color_theme_name: 'element-symbol',
            molstar_color_theme_params: { carbonColor: options?.carbon_color ? { name: 'uniform', params: { value: decodeColor(options?.carbon_color) } } : { name: 'element-symbol', params: { } } }
        }
    };
    comp.representation({ type: 'surface' }).color(coloring).opacity({ opacity: 0.33 });

    return comp;
}

function bindingSite(structure: MVSStructure, residues: [selector: ComponentExpressionT, label: string][], options: {
    color?: ColorT,
    label_size?: number,
    label_color?: ColorT,
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
            label_color: options?.label_color ?? color,
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
            linger_duration_ms: 5000,
            transition_duration_ms: 1500,
        });
    });

    return {
        kind: 'multiple',
        snapshots,
        metadata: {
            title: 'The Structural Story of TATA-Binding Protein and its Role in Transcription Initiation',
            version: '1.0',
            timestamp: new Date().toISOString(),
        }
    };
}