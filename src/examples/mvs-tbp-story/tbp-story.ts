/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { MVSData_States } from '../../extensions/mvs/mvs-data';
import { createMVSBuilder, Structure as MVSStructure, Root } from '../../extensions/mvs/tree/mvs/mvs-builder';
import { MVSNodeParams } from '../../extensions/mvs/tree/mvs/mvs-tree';
import { ColorT } from '../../extensions/mvs/tree/mvs/param-types';
import { Mat3, Mat4, Vec3 } from '../../mol-math/linear-algebra';

const Colors = {
    '1vok': '#4577B2' as ColorT,
    '1cdw': '#BC536D' as ColorT,
    '1vtl': '#B9E3A0' as ColorT,
    '7enc': '#F3774B' as ColorT,

    'active-site': '#F3794C' as ColorT,
    'binding-site': '#FEEB9F' as ColorT,
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
The TATA-binding protein (TBP) is the central element of this system.

The crystal structure of TBP from *Arabidopsis thaliana* in its apo form is shown ([PDB ID 1VOK](${wwPDBLink('1vok')})).
The structure reveals a highly symmetric DNA-binding fold composed of two topologically identical domains derived from the two direct repeats in the phylogenetically conserved sequence.
The domains are related by an approximate 2-fold axis, and consist of a five-stranded, antiparallel beta sheet and two alpha helices.
It is thought that an ancient gene duplication created this protein by combining two copies of the same gene. 

The intramolecular symmetry generates a saddle-shaped structure dominated by a curved antiparallel beta sheet, which forms the saddle's concave face and interacts with the DNA.
The convex face of the saddle would interact with other proteins during transcription.
`,
        state: (): Root => {
            const builder = createMVSBuilder();

            const _1vok = structure(builder, '1vok');
            const [_1vok_poly,] = select(_1vok, { color: Colors['1vok'] });
            _1vok_poly.label({ text: 'TATA-Binding Protein' });

            return builder;
        },
        camera: {
            position: [155.18, 118.49, 49.18],
            target: [74.81, 66.8, 30.7],
            up: [-0.55, 0.83, 0.1],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'TBP: Highly Conserved in Eukaryotes',
        key: 'highly-conserved',
        description: `
### TATA-Binding Protein Is Highly Conserved in Eukaryotes

Eukaryotic protein-coding genes have a characteristic sequence of nucleotides, termed the TATA box, in front of the start site of transcription.
The typical sequence is something like T-A-T-A-a/t-A-a/t, where a/t refers to positions that can be either A or T.
TBP recognizes this TATA sequence and binds to it, creating a landmark that marks the start site of transcription. 

TBP has a phylogenetically conserved, 180 amino-acid carboxy-terminal portion (ranging from 38 to 93% identity), containing two structural repeats flanking a highly basic segment known as the basic repeat.
The C-terminal or core portion of the protein binds the TATA consensus sequence with high affinity, interacting with the minor groove and promoting DNA bending.
Structural superposition of TBP bound to DNA in human ([PDB ID 1CDW](${wwPDBLink('1cdw')})) and *Arabidopsis thaliana* ([PDB ID 1VTL](${wwPDBLink('1vtl')})) shows that their sequences are 83% identical and the RMSD between the structures is 0.43 Å.
The structures of *Arabidopsis thaliana* and human core TBP-TATA element co-crystal structures demonstrate a common induced-fit mechanism of protein-DNA recognition involving subtle conformation changes in the protein and an unprecedented DNA distortion. 
`,
        state: () => {
            const builder = createMVSBuilder();

            const _1cdw = structure(builder, '1cdw');
            select(_1cdw, { color: Colors['1cdw'], selector: 'protein' });
            select(_1cdw, { color: Colors['1cdw'], selector: 'nucleic', opacity: 0.5 });
            // TODO domain for TATA fragment?

            const _1vtl = structure(builder, '1vtl');
            select(_1vtl, { color: Colors['1vtl'], selector: { label_asym_id: 'E' } });
            select(_1vtl, { color: Colors['1vtl'], selector: { label_asym_id: 'A' }, opacity: 0.5 });
            select(_1vtl, { color: Colors['1vtl'], selector: { label_asym_id: 'B' }, opacity: 0.5 });
            // TODO domain for TATA fragment?

            return builder;
        },
        camera: {
            position: [155.18, 118.49, 49.18],
            target: [74.81, 66.8, 30.7],
            up: [-0.55, 0.83, 0.1],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'TBP: Binding to TATA Box [1/5]',
        key: 'tata-box-overview',
        transition_duration_ms: 750,
        description: `
### TATA-Binding Protein Binds to the Minor Groove of the TATA Box

When the first structures of TBP were determined, researchers discovered that TBP is not gentle when it binds to DNA.
Instead, it grabs the TATA sequence and bends it sharply, as seen in PDB entries ([PDB ID 1VTL](${wwPDBLink('1vtl')})) and ([PDB ID 1CDW](${wwPDBLink('1cdw')})). 

Interactions with the minor groove can be divided into different classes as seen in PDB ([PDB ID 1CDW](${wwPDBLink('1cdw')})).
`,
        state: () => {
            const builder = createMVSBuilder();

            const _1cdw = structure(builder, '1cdw');
            select(_1cdw, { color: Colors['1cdw'], selector: { label_asym_id: 'C' } });
            select(_1cdw, { color: Colors['1cdw'], selector: { label_asym_id: 'A' } });
            select(_1cdw, { color: Colors['1cdw'], selector: { label_asym_id: 'B' } });

            return builder;
        },
        camera: {
            position: [155.18, 118.49, 49.18],
            target: [74.81, 66.8, 30.7],
            up: [-0.55, 0.83, 0.1],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'TBP: Basic Amino Acids [2/5]',
        key: 'tata-box-1',
        transition_duration_ms: 750,
        description: `
### Basic Amino Acids

A string of lysine and arginine amino acids interact with the phosphate groups of the DNA (either directly or mediated through water) and glues the protein to the DNA (Arg192, Arg199, Arg290, Lys221, Lys312, Lys204, Lys214, Lys295, Lys305).
`,
        state: () => {
            const builder = createMVSBuilder();

            const _1cdw = structure(builder, '1cdw');
            select(_1cdw, { color: Colors['1cdw'], selector: { label_asym_id: 'C' } });
            select(_1cdw, { color: Colors['1cdw'], selector: { label_asym_id: 'A' } });
            select(_1cdw, { color: Colors['1cdw'], selector: { label_asym_id: 'B' } });
            // Arg: 38, 45, 136
            // Lys: 67, 158, 50, 60, 141, 151

            return builder;
        },
        camera: {
            position: [155.18, 118.49, 49.18],
            target: [74.81, 66.8, 30.7],
            up: [-0.55, 0.83, 0.1],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'TBP: Phenylalanine [3/5]',
        key: 'tata-box-2',
        transition_duration_ms: 750,
        description: `
### Phenylalanine-Induced Kinks

At either end of the TATA element there are two pairs of phenylalanine side chains partially inserted between adjacent base pairs, producing the two dramatic kinks (Phe193, Phe210, Phe284, and Phe301). Between the two kinks the DNA is partially unwound.
`,
        state: () => {
            const builder = createMVSBuilder();

            const _1cdw = structure(builder, '1cdw');
            select(_1cdw, { color: Colors['1cdw'], selector: { label_asym_id: 'C' } });
            select(_1cdw, { color: Colors['1cdw'], selector: { label_asym_id: 'A' } });
            select(_1cdw, { color: Colors['1cdw'], selector: { label_asym_id: 'B' } });
            // Phe: 39, 56, 130, 147

            return builder;
        },
        camera: {
            position: [155.18, 118.49, 49.18],
            target: [74.81, 66.8, 30.7],
            up: [-0.55, 0.83, 0.1],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'TBP: H-Bonds in Minor Groove [4/5]',
        key: 'tata-box-3',
        transition_duration_ms: 750,
        description: `
### Hydrogen Bonds in the Minor Groove

Polar side chains make minor groove hydrogen bonds with acceptors of base pairs centred about the approximate 2-fold symmetry axis (Asn163, Asn253, and Thr309).
`,
        state: () => {
            const builder = createMVSBuilder();

            const _1cdw = structure(builder, '1cdw');
            select(_1cdw, { color: Colors['1cdw'], selector: { label_asym_id: 'C' } });
            select(_1cdw, { color: Colors['1cdw'], selector: { label_asym_id: 'A' } });
            select(_1cdw, { color: Colors['1cdw'], selector: { label_asym_id: 'B' } });
            // Asn: 9, 99, 155

            return builder;
        },
        camera: {
            position: [155.18, 118.49, 49.18],
            target: [74.81, 66.8, 30.7],
            up: [-0.55, 0.83, 0.1],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'TBP: Non-Polar Interactions [5/5]',
        key: 'tata-box-4',
        transition_duration_ms: 750,
        description: `
### Hydrophobic and van der Waals Interactions

Several residues projecting from the concave surface of TBP make hydrophobic or van der Waals side-chain/base contacts (<4 Å between non-hydrogen atoms) with the TATA element. The combination of interactions allows TATA-binding protein to recognize the proper DNA sequence.
`,
        state: () => {
            const builder = createMVSBuilder();

            const _1cdw = structure(builder, '1cdw');
            select(_1cdw, { color: Colors['1cdw'], selector: { label_asym_id: 'C' } });
            select(_1cdw, { color: Colors['1cdw'], selector: { label_asym_id: 'A' } });
            select(_1cdw, { color: Colors['1cdw'], selector: { label_asym_id: 'B' } });
            // TODO need selection

            return builder;
        },
        camera: {
            position: [155.18, 118.49, 49.18],
            target: [74.81, 66.8, 30.7],
            up: [-0.55, 0.83, 0.1],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'TBP and Transcription Pre-Initiation Complex',
        key: 'pre-init-complex',
        description: `
### TATA-Binding Protein and the Transcription Pre-Initiation Complex

The structure shown here, from [PDB ID 7ENC](${wwPDBLink('7enc')}), is a “pre-initiation complex” (PIC) poised to start transcription.
This complex includes eukaryotic Mediator and RNA polymerase II, along with a collection of general transcription factors that perform the necessary tasks of recognizing the transcription start site of genes, separating the two strands of the DNA double helix, and facilitating transcription initiation by RNA polymerase II.

TBP starts the process of transcription and works as part of a larger transcription factor, TFIID, that belongs to the collection of general transcription factors that cradle RNA polymerase in the PIC.
After TBP binds to the promoter, it recruits additional transcription factors.
TFIIA and TFIIB interact with surrounding regions of the DNA and, along with TFIIF, assist with positioning of RNA polymerase II at the transcription start site.
TFIIE, TFIIH, and other TFIID subunits bring additional functionality to the complex.
In particular, part of TFIIH is a translocase that separates the two strands of DNA in preparation for transcription, and the CAK module of TFIIH adds phosphate groups to the long tail of RNA polymerase II, sending the signal that it is time to get started with mRNA synthesis. 

TBP is a critical player that gets the transcription process started and is central to the pre-initiation complex.
Under tremendous selection pressure, it is highly conserved in eukaryotes.
`,
        state: () => {
            const builder = createMVSBuilder();

            const _7enc = structure(builder, '7enc');
            select(_7enc, { color: Colors['7enc'], selector: { label_asym_id: 'CB' } });
            select(_7enc, { color: Colors['7enc'], selector: { label_asym_id: 'GB' } });
            select(_7enc, { color: Colors['7enc'], selector: { label_asym_id: 'HB' } });
            select(_7enc, { color: Colors['7enc'], opacity: 0.5 });

            return builder;
        },
        camera: {
            position: [155.18, 118.49, 49.18],
            target: [74.81, 66.8, 30.7],
            up: [-0.55, 0.83, 0.1],
        } satisfies MVSNodeParams<'camera'>,
    }, {
        header: 'The End',
        key: 'end',
        description: `
### The End

That's all folks! We hope you enjoyed this interactive journey through the structural biology of the TATA-binding protein.

The next time you look at a macromolecular structure, remember: each atom tells a story, and each discovery shapes the future of medicine.

Read more [here](https://pdb101.rcsb.org/motm/289).
`,
        state: (): Root => {
            return Steps[Steps.length - 2].state();
        },
        camera: {
            position: [591.05, 316.22, 162.66],
            target: [82.39, -10.93, 45.7],
            up: [-0.55, 0.83, 0.1],
        } satisfies MVSNodeParams<'camera'>,
    }
];

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
            transition_duration_ms: s.transition_duration_ms ?? 1500,
        });
    });

    return {
        kind: 'multiple',
        snapshots,
        metadata: {
            title: 'The Structural Story of TATA-Binding Protein and Its Role in Transcription Initiation',
            version: '1.0',
            timestamp: new Date().toISOString(),
        }
    };
}