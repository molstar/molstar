/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { decodeColor } from '../../../extensions/mvs/helpers/utils';
import { MVSData_States } from '../../../extensions/mvs/mvs-data';
import { createMVSBuilder, Structure as MVSStructure, Root } from '../../../extensions/mvs/tree/mvs/mvs-builder';
import { MVSNodeParams } from '../../../extensions/mvs/tree/mvs/mvs-tree';
import { ColorT } from '../../../extensions/mvs/tree/mvs/param-types';
import { Mat4 } from '../../../mol-math/linear-algebra/3d/mat4';
import { Vec3 } from '../../../mol-math/linear-algebra/3d/vec3';

// 1pmb->1mbn
const align = Mat4.fromArray(Mat4.zero(), [0.4634187130865737,-0.7131589697034304,0.5259728687171936,0,-0.22944227902330105,-0.6698811108214233,-0.7061273127008398,0,0.8559202154942049,0.2065522332899299,-0.4740643150728161,0,-52.54880970106205,37.49099778180445,-6.133850309914719,1], 0);
const translation = Mat4.fromTranslation(Mat4.zero(), Vec3.create(50, 0, 0));

const Colors = {
    '1mbn': '#b26445' as ColorT,

    'ligand-away': '#F3794C' as ColorT,
    'ligand-docked': '#B9E3A0' as ColorT,
};

const GColors = {
    molstar_color_theme_name: 'element-symbol',
    molstar_color_theme_params: {
        carbonColor: {
            name: 'uniform',
            params: { value: decodeColor('#ffffff') }
        },
    }
};

const GColors2 = {
                    molstar_color_theme_name: 'illustrative',
                    molstar_color_theme_params: {
                        carbonColor: {
                            name: 'uniform',
                            params: { value: decodeColor('#bd9e9eff') }
                        },
                    }
                };
/* from David Goodsell style
in his illustrate software
HETATM-H-------- 0,9999, 1.1,1.1,1.1, 0.0
HETATMH--------- 0,9999, 1.0,1.0,1.0, 0.0
ATOM  -H-------- 0,9999, 1.0,1.0,1.0, 0.0
ATOM  H--------- 0,9999, 1.0,1.0,1.0, 0.0
HETATM-----HOH-- 0,9999, 1.0,1.0,0.0, 0.0
ATOM  -OD--ASP A 0,9999  1.00, 0.20, 0.20, 1.6
ATOM  -OE--GLU A 0,9999  1.00, 0.20, 0.20, 1.6
ATOM  -NZ--LYS A 0,9999  0.10, 0.70, 1.00, 1.6
ATOM  -NH--ARG A 0,9999  0.10, 0.70, 1.00, 1.6
ATOM  -NE--ARG A 0,9999  0.10, 0.70, 1.00, 1.6
ATOM  -ND--HIS A 0,9999  0.10, 0.70, 1.00, 1.6
ATOM  -NE--HIS A 0,9999  0.10, 0.70, 1.00, 1.6
ATOM  -N------ A 0,9999  0.80, 0.90, 1.00, 1.5
ATOM  -O------ A 0,9999  1.00, 0.80, 0.80, 1.5
ATOM  -C------ A 0,9999  1.00, 1.00, 1.00, 1.6
ATOM  -S------ A 0,9999  1.00, 0.90, 0.50, 1.8
HETATM-C------ - 0,9999  0.60, 0.90, 0.60, 1.5
HETATM-------- - 0,9999  0.40, 0.90, 0.40, 1.5
*/
const GColors3 = {
                schema: 'all_atomic', // or maybe just 'atom'
                category_name: 'atom_site',
                field_name: 'type_symbol',
                palette: {
                    kind: 'categorical',
                    // missing_color: ...
                    colors: {
                        'C': '#FFFFFF',
                        'N': '#CCE6FF',
                        'O': '#FFCCCC',
                        'S': '#FFE680',
                    }
                }
            } as unknown as MVSNodeParams<'color_from_source'>;
// Python script
// import base64
// def convert(fname):
//     with open(fname, "rb") as f:
//         b64 = base64.b64encode(f.read()).decode("utf-8")
//     print(f"data:audio/mp4;base64,{b64}")
//     return b64 


const _Audio1 = "https://raw.githubusercontent.com/molstar/molstar/master/examples/audio/AudioMOM1_A.mp3";
const _Audio2 = "https://raw.githubusercontent.com/molstar/molstar/master/examples/audio/AudioMOM1_B.mp3";
const _Audio3 = "https://raw.githubusercontent.com/molstar/molstar/master/examples/audio/AudioMOM1_C.mp3";
const _Audio4 = "https://raw.githubusercontent.com/molstar/molstar/master/examples/audio/AudioMOM1_D.mp3";

const description_intro = `
# Molecule of the Month: Myoglobin
Basic controls for the audio comments :
[Play](${encodeURIComponent(`!play-audio=${_Audio1}`)})
[Pause](!pause-audio)
[Stop](!stop-audio)

Myoglobin was the first protein to have its atomic structure determined, revealing how it stores oxygen in muscle cells.

---

## The First Protein Structure

Any discussion of protein structure must necessarily begin with **myoglobin**, because it is where the science of protein structure began. After years of arduous work, *John Kendrew* and his coworkers determined the atomic structure of myoglobin, laying the foundation for an era of biological understanding.

You can take a close look at this protein structure yourself, in **PDB entry [1mbn](https://www.rcsb.org/structure/1MBN)**. You will be amazed—just like the world was in 1960—at the beautiful intricacy of this protein.

---

## Myoglobin and Muscles

[Myoglobin](!query%3Dchain%20A%26lang%3Dpymol%26action%3Dhighlight%2Cfocus) is a **small, bright red protein**. It is very common in muscle cells and gives meat much of its red color. Its job is to **store oxygen**, for use when muscles are hard at work.

To do this, it uses a special chemical tool to capture slippery oxygen molecules: a **[heme group](!query%3Dresn%20HEM%26lang%3Dpymol%26action%3Dhighlight%2Cfocus)**. Heme is a disk-shaped molecule with a hole in the center that is perfect for holding an iron ion. The iron then forms a strong interaction with the **[oxygen molecule](!query%3Dresn%20OH%26lang%3Dpymol%26action%3Dhighlight%2Cfocus)**. As you can see in the structure, the heme group is held tightly in a deep pocket on one side of the protein.

---

## Visualizing Protein Structure

When the structure of myoglobin was solved, it posed a great challenge. The structure is so complex that **new methods** needed to be developed to display and understand it.

- *John Kendrew* used a huge wire model to build the structure based on the experimental electron density.  
- Then, the artist *Irving Geis* was employed to create a picture of myoglobin for a prominent article in *Scientific American*.  
- Computer graphics were still many years in the future, so he created this illustration entirely by hand—one atom at a time.  

You can learn more about the work of Irving Geis at the **[Geis Archive on PDB-101](https://pdb101.rcsb.org/learn/GeisArchive)**.

![Alt Text](https://cdn.rcsb.org/pdb101/motm/1/1-Myoglobin-geis-0218-myoglobin.png)
*Illustration of myoglobin by Irving Geis. You can learn more about this painting at the Geis Archive on PDB-101.  
Used with permission from the Howard Hughes Medical Institute, Copyright 2015.*


`;

const description_p1=`
# Myoglobin and Whales
Basic controls for the audio comments :
[Play](${encodeURIComponent(`!play-audio=${_Audio2}`)})
[Pause](!pause-audio)
[Stop](!stop-audio)

If you look at John Kendrew's PDB file, you'll notice that the myoglobin he used was taken 
from sperm whale muscles. Whales and dolphin have a great need for myoglobin, so that they can 
store extra oxygen for use in their deep undersea dives. Typically, they have about 30 
times more than in animals that live on land. A recent study revealed that a few special 
modifications are needed to make this possible. Comparing [whale myoglobin](!query%3Dmodel%201mbn%26lang%3Dpymol%26action%3Dhighlight%2Cfocus) (PDB entry [1mbn](https://www.rcsb.org/structure/1mbn)) with 
[pig myoglobin](!query%3Dmodel%201pmb%26lang%3Dpymol%26action%3Dhighlight%2Cfocus) (PDB entry [1pmb](https://www.rcsb.org/structure/1pmb)), we find that there are several mutations that add extra positively-charged 
amino acids to the surface. Marine animals typically have these extra charges on the surface of their myoglobin 
to help repel neighboring molecules and prevent aggregation when myoglobin is at high concentrations.
`;

const description_p2=`
# Oxygen Bound to Myoglobin
Basic controls for the audio comments :
[Play](${encodeURIComponent(`!play-audio=${_Audio3}`)})
[Pause](!pause-audio)
[Stop](!stop-audio)

A later structure of myoglobin, PDB entry 1mbo, shows that oxygen binds to 
the iron atom deep inside the protein. So how does it get in and out? The 
answer is that the structure in the PDB is only one snapshot of the 
protein, caught when it is in a tightly-closed form. In reality, 
myoglobin (and all other proteins) is constantly in motion, performing 
small flexing and breathing motions. So, temporary openings constantly 
appear and disappear, allowing oxygen in and out.
`;

const description_p3 = `
# Molecule of the Month: Myoglobin
Basic controls for the audio comments :
[Play](${encodeURIComponent(`!play-audio=${_Audio1}`)})
[Pause](!pause-audio)
[Stop](!stop-audio)

The atomic structure of myoglobin revealed many of the basic principles 
of protein structure and stability. For instance, the structure showed 
that when the protein chain folds into a globular structure, 
[carbon-rich amino acids]() are sheltered inside and [charged amino acids]() 
are most often found on the surface, occasionally forming salt bridges 
that pair two opposite charges (shown here with circles). 

To explore some of these principles, eplxore freely in the interactive view.

# Topics for Further Discussion
You can use the sequence comparison tool to align the sequences of different myoglobins, looking for mutations. For instance, here is the alignment of whale and pig myoglobin used to create the illustration in this column.
PDB entry 2jho includes myoglobin poisoned by cyanide. Take a look and you'll see that the cyanide blocks the binding site for oxygen.


# References

- 1mbn: J. C. Kendrew, R. E. Dickerson, B. E. Strandberg, R. G. Hart, D. R. Davies, D. C. Phillips & V. C. Shore (1960) Structure of Myoglobin. Nature 185, 422-427.
- J. C. Kendrew (1961) The three-dimensional structure of a protein molecule. Scientific American 205(6), 96-110.
- 1mbo: S. E. Phillips (1980) Structure and refinement of oxymyoglobin at 1.6 A resolution. Journal of Molecular Biology 142, 531-554.
- 1pmb: S. J. Smerdon, T. J. Oldfield, E. J. Dodson, G. G. Dodson, R. E. Hubbard & A. J. Wilkinson (1990) Determination of the crystal structure of recombinant pig myoglobin by molecular replacement and its refinement. Acta Crystallographica B, 46, 370-377.
- S. Mirceta, A. V. Signore, J. M. Burns, A. R. Cossins, K. L. Campbell & M. Berenbrink (2013) Evolution of mammalian diving capacity traced by myoglobin net surface charge. Science 340, 1234192.

`;

const Steps = [
    {
        header: 'Molecule of the Month: Myoglobin',
        key: 'intro',
        description: description_intro,
        linger_duration_ms: 2000,
        transition_duration_ms: 500,
        state: (): Root => {
            const builder = createMVSBuilder();

            const _1mbn = structure(builder, '1MBN');
            polymer(_1mbn, { color: Colors['1mbn'] });
            _1mbn.component({ selector: 'ligand' })
            .representation({ ref: 'ligand', type: 'ball_and_stick',
                custom: {
                    molstar_reprepresentation_params: {
                        emissive: 0.0
                    }
                }
            })
            .color({ color: 'orange' });
            // FE and O should be spacefill

            _1mbn.component({ selector: { auth_seq_id: 155, label_atom_id: 'F' } })
            .representation({ type: 'spacefill' })
            .color({ color: 'yellow' });

            _1mbn.component({ selector: { auth_seq_id: 154 } })
            .representation({ type: 'spacefill' })
            .color({ color: 'blue' });

            const chA = _1mbn.component({ selector: { label_asym_id: 'A' } });
            chA.representation({ type: 'surface', surface_type: 'gaussian' })
            .color({ color: '#ff0303' })
            .opacity({ ref: 'surfopa', opacity: 0.0 });

            chA.representation({ type: 'line' })
            .color({ custom: { molstar_color_theme_name: "element-symbol" } })
            .opacity({ ref: 'lineopa', opacity: 0.0 });

            const prims = _1mbn.primitives({
                ref: 'prims',
                label_opacity: 1,
                label_background_color: 'grey',
                custom: {
                    molstar_markdown_commands: {
                        // 'apply-snapshot': 'interlude',
                        'play-audio': _Audio1,
                    }
                }
            });
            prims.label({
                text: 'Start Comments',
                position: [13.5, 45.1, 7.7],
                label_size: 10
            });
            // doesnt work for first slide, but work afterward
            builder.extendRootCustomState({
                molstar_on_load_markdown_commands: {
                    'play-audio': _Audio1,
                }
            });
            const anim = builder.animation(
                {
                    custom: {
                    molstar_trackball: {
                        name: 'spin',
                        params: { speed: -0.05 },
                    }
                }
            });
            anim.interpolate({
                kind: 'scalar',
                target_ref: 'lineopa',
                duration_ms: 2000,
                start_ms: 0,
                // frequency: 4,
                // alternate_direction: true,
                property: 'opacity',
                start: 0.0,
                end: 1.0,
            });
            // anim.interpolate({
            //     kind: 'scalar',
            //     target_ref: 'surfopa',
            //     duration_ms: 2000,
            //     start_ms: 0,
            //     // frequency: 4,
            //     // alternate_direction: true,
            //     property: 'opacity',
            //     start: 1.0,
            //     end: 0.0,
            // });
            // // emission of hem group ligand
            // anim.interpolate({
            //     kind: 'scalar',
            //     target_ref: 'ligand',
            //     start_ms: 8000,
            //     duration_ms: 2000,
            //     frequency: 4,
            //     alternate_direction: true,
            //     property: ['custom', 'molstar_reprepresentation_params', 'emissive'],
            //     end: 1.0,
            // });
            return builder;
        },
        camera: {
            position: [13.5, 21.1, 73.1],
            target: [13.5, 21.1, 7.7],
            up: [0,1,0],
        } satisfies MVSNodeParams<'camera'>,
    },
    {
        header: 'Myoglobin and Whales',
        key: 'whale',
        description: description_p1,
        linger_duration_ms: 2000,
        transition_duration_ms: 500,
        state: (): Root => {
            const builder = createMVSBuilder();
            // use alignement and translation ?
            const _1mbn = structure(builder, '1mbn').transform({ translation: [-60, 0, 0] });
            // TODO: align the two structure and animate the separation
            // to show the difference in charged amino acids
            // whale
            _1mbn.component({ selector: { label_asym_id: 'A' } })
            .representation({ type: 'spacefill', custom: { molstar_representation_params: { ignoreLight: true } } })
            .colorFromSource({
                schema: 'all_atomic', // or maybe just 'atom'
                category_name: 'atom_site',
                field_name: 'type_symbol',
                palette: {
                    kind: 'categorical',
                    // missing_color: ...
                    colors: {
                        'C': '#FFFFFF',
                        'N': '#CCE6FF',
                        'O': '#FFCCCC',
                        'S': '#FFE680',
                    }
                }
            });

            _1mbn.component({ selector: { auth_seq_id: 155 } })
            .representation({ type: 'spacefill', custom: { molstar_representation_params: { ignoreLight: true } } })
            .color({ custom: GColors2 });

            _1mbn.primitives({
                ref: 'prims',
                label_opacity: 1,
                label_attachment: 'top-center',
                label_show_tether: true,
                label_tether_length: 3.0,
                custom: {
                    molstar_markdown_commands: {
                        // 'stop-audio': true,
                        // what can we do here ?
                    }
                }
            })
            .label({ text: 'whale',
                    position: { label_asym_id: 'A' },
                    label_size: 10 });

            _1mbn.primitives({
                label_opacity: 1,
                // label_color: '#000000',
             })
            .label({ text: '★',
                label_offset: 4,
                position: { label_asym_id: 'A', auth_seq_id: 12, atom_id: 96 }, label_size: 5 })
            .label({ text: '★',label_offset: 4,
                position: { label_asym_id: 'A', auth_seq_id: 140, auth_atom_id: 'NZ' }, label_size: 5 })
            .label({ text: '★',label_offset: 4,
                position: { label_asym_id: 'A', auth_seq_id: 87, auth_atom_id: 'NZ' }, label_size: 5 });

            const _1pmb = structure(builder, '1pmb').transform({ matrix: align });

            _1pmb.component({ selector: { label_asym_id: 'A' } })
            .representation({ type: 'spacefill' , custom: { molstar_representation_params: { ignoreLight: true } } })
            .colorFromSource(GColors3);

            _1pmb.component({ selector: { label_asym_id: 'C', auth_seq_id: 154 } })
            .representation({ type: 'spacefill' , custom: { molstar_representation_params: { ignoreLight: true } } })
            .color({ custom: GColors2 });

            _1pmb.primitives({
                ref: 'prims',
                label_opacity: 1,
                label_attachment: 'top-center',
                label_show_tether: true,
                label_tether_length: 3.0,
                custom: {
                    molstar_markdown_commands: {
                        // 'stop-audio': true,
                        // what can we do here ?
                    }
                }
            })
            .label({ text: 'pig',
                    position: { label_asym_id: 'A' },
                    label_size: 10 });

            builder.extendRootCustomState({
                molstar_on_load_markdown_commands: {
                    'play-audio': _Audio2,
                }
            });

            // const prims = _1mbn.primitives({
            //     ref: 'prims',
            //     label_opacity: 1,
            //     custom: {
            //         molstar_markdown_commands: {
            //             'stop-audio': true,
            //         }
            //     }
            // });
            // prims.label({ text: 'Stop', position: { label_asym_id: 'A' }, label_size: 10 });
            const anim = builder.animation(
                {
                    custom: {
                    molstar_trackball: {
                        name: 'spin',
                        params: { speed: -0.05 },
                    }
                }
            });
            return builder;
        },
        camera: {
            position: [-14.6, 116.1, 66.5],
            target: [-18.9, 21.1, 7.3],
            up: [-0.0, 0.5, -0.8],
        } satisfies MVSNodeParams<'camera'>,
    },
     {
        header: 'Oxygen Bound',
        key: 'oxygen',
        description: description_p2,
        linger_duration_ms: 2000,
        transition_duration_ms: 500,
        state: (): Root => {
            const builder = createMVSBuilder();
            const _1mbo = structure(builder, '1mbo');
            // TODO: find a way to illustrate the breathing
            _1mbo.component({ selector: { label_asym_id: 'A' } })
            .transform({ translation: [0, 0, 0] })
            .representation({ type: 'spacefill' })
            // .color({ custom: { molstar_color_theme_name: "illustrative" } })
            .color({ custom: GColors2 })
            .opacity({ ref: 'spo', opacity: 1.0 });

            _1mbo.component({ selector: { label_asym_id: 'C', auth_seq_id: 155 } })
            .representation({ type: 'spacefill' })
            .color({ custom: GColors2 });

            _1mbo.component({ selector: { label_asym_id: 'A' } })
            .representation({ type: 'backbone' })
            .color({ color: '#919191' });

            _1mbo.component({ selector: { label_asym_id: 'D', auth_seq_id: 555 } })
            .representation({ type: 'spacefill' })
            .color({ color: '#8670f0' });

            builder.extendRootCustomState({
                molstar_on_load_markdown_commands: {
                    'play-audio': _Audio3,
                }
            });
            const anim = builder.animation(
                {
                    custom: {
                    molstar_trackball: {
                        name: 'spin',
                        params: { speed: -0.05 },
                    }
                }
            });
            anim.interpolate({
                kind: 'scalar',
                target_ref: 'spo',
                duration_ms: 5000,
                start_ms: 0,
                // frequency: 4,
                // alternate_direction: true,
                property: 'opacity',
                start: 1.0,
                end: 0.05,
            });
            return builder;
        },
        camera: {
            position: [13.5, 21.1, 73.1],
            target: [13.5, 21.1, 7.7],
            up: [0,1,0],
        } satisfies MVSNodeParams<'camera'>,
    },
     {
        header: 'Conclusion',
        key: 'end',
        description: description_p3,
        linger_duration_ms: 2000,
        transition_duration_ms: 500,
        state: (): Root => {
            const builder = createMVSBuilder();
            const _1mbn = structure(builder, '1mbn');
            // TODO:
            // select all carbon rich
            // resn ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO
            const carb = ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO'].map(amk => ({ label_comp_id: amk }));
            // resn LYS+ARG+HIS+ASP+GLU
            const chargedp = ['LYS', 'ARG', 'HIS'].map(amk => ({ label_comp_id: amk }));
            const chargedn = ['ASP', 'GLU'].map(amk => ({ label_comp_id: amk }));
            // select all charged residues
            const a = _1mbn.component({ selector: carb });
            a.representation({ type: 'ball_and_stick' })
            .color({ color: '#c8c8c7' });

            const b = _1mbn.component({ selector: chargedp });
            b.representation({ type: 'ball_and_stick' })
            .color({ color: '#ca2f2f' });

            const c = _1mbn.component({ selector: chargedn });
            c.representation({ type: 'ball_and_stick' })
            .color({ color: '#6273ce' });

            _1mbn.component({ selector: { label_asym_id: 'A' } })
            .representation({ type: 'backbone' })
            .color({ color: '#919191' });

            _1mbn.component({ selector: 'ligand' })
            .representation({ ref: 'ligand', type: 'ball_and_stick',
                custom: {
                    molstar_reprepresentation_params: {
                        emissive: 0.0
                    }
                }
            })
            .color({ color: 'orange' });

            builder.extendRootCustomState({
                molstar_on_load_markdown_commands: {
                    'play-audio': _Audio4,
                }
            });

            const anim = builder.animation({});
            return builder;
        },
        camera: {
            position: [16.0, 47.2, 67.8],
            target: [13.6, 21.1, 7.6],
            up: [0.1, 0.9, -0.4],
        } satisfies MVSNodeParams<'camera'>,
    },
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
            title: 'Audio Showcase',
            version: '1.0',
            timestamp: new Date().toISOString(),
        }
    };
}