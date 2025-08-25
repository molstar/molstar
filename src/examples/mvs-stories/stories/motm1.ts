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
import { Mat4 } from '../../../mol-math/linear-algebra/3d/mat4';
import { Vec3 } from '../../../mol-math/linear-algebra/3d/vec3';
import { MolScriptBuilder as MS } from '../../../mol-script/language/builder';
import { formatMolScript } from '../../../mol-script/language/expression-formatter';

// 1pmb->1mbn
const align = Mat4.fromArray(Mat4.zero(), [0.4634187130865737, -0.7131589697034304, 0.5259728687171936, 0, -0.22944227902330105, -0.6698811108214233, -0.7061273127008398, 0, 0.8559202154942049, 0.2065522332899299, -0.4740643150728161, 0, -52.54880970106205, 37.49099778180445, -6.133850309914719, 1], 0);
// 1mbo->1myf
const alignmbo = Mat4.fromArray(Mat4.zero(), [-0.8334619943964441, -0.512838061396133, -0.20576353166796402, 0, -0.20145089001561267, 0.628743285359846, -0.7510655776229758, 0, 0.5145474196737698, -0.5845332204089626, -0.6273453801378679, 0, 11.864847328611186, -1.5261713438028912, 23.638919347623467, 1], 0);

const ill_color = (color: string, carbonLightness: number) => ({
    molstar_color_theme_name: 'illustrative',
    molstar_color_theme_params: {
        style: {
            name: 'uniform',
            params: {
                value: decodeColor(color),
                saturation: 0,
                lightness: 0,
            }
        },
        carbonLightness: carbonLightness // required parameter
    }
});

const GColors2 = ill_color('#947c7c', 0.8);

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

// const path = "https://raw.githubusercontent.com/molstar/molstar/master";
const path = '';
const _Audio1 = path + '/examples/audio/AudioMOM1_A.mp3';
const _Audio2 = path + '/examples/audio/AudioMOM1_B.mp3';
const _Audio3 = path + '/examples/audio/AudioMOM1_C.mp3';
const _Audio4 = path + '/examples/audio/AudioMOM1_D.mp3';

const q = (expr: string, lang = 'pymol') =>
    `!query=${encodeURIComponent(expr)}&lang=${lang}&action=highlight,focus`;

const description_intro = `
# Molecule of the Month: Myoglobin
A story based on the orginal [first Molecule of the Month](https://pdb101.rcsb.org/motm/1) made by David Goodsell in January 2000.


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


const query1 = MS.struct.generator.atomGroups({
    'entity-test': MS.core.rel.eq([
        MS.struct.atomProperty.core.modelEntryId(),
        '1MBN'
    ])
});
const firstEntity1 = q(formatMolScript(query1), 'mol-script');
const query2 = MS.struct.generator.atomGroups({
    'entity-test': MS.core.rel.eq([
        MS.struct.atomProperty.core.modelEntryId(),
        '1PMB'
    ])
});
const firstEntity2 = q(formatMolScript(query2), 'mol-script');

const query3 = MS.struct.generator.atomGroups({
    'entity-test': MS.core.rel.eq([
        MS.struct.atomProperty.core.modelEntryId(),
        '1MBN'
    ]),
    'residue-test': MS.core.set.has([
        MS.set(12, 140, 87),
        MS.struct.atomProperty.macromolecular.auth_seq_id()
    ])
});
const charged_residues = q(formatMolScript(query3), 'mol-script');

const description_p1 = `
# Myoglobin and Whales
Basic controls for the audio comments :
[Play](${encodeURIComponent(`!play-audio=${_Audio2}`)})
[Pause](!pause-audio)
[Stop](!stop-audio)

If you look at John Kendrew's PDB file, you'll notice that the myoglobin he used was taken 
from sperm whale muscles. Whales and dolphin have a great need for myoglobin, so that they can 
store extra oxygen for use in their deep undersea dives. Typically, they have about 30 
times more than in animals that live on land. A recent study revealed that a few special 
modifications are needed to make this possible.
Comparing [whale myoglobin](${firstEntity1})
(PDB entry [1mbn](https://www.rcsb.org/structure/1mbn)) with 
[pig myoglobin](${firstEntity2})
(PDB entry [1pmb](https://www.rcsb.org/structure/1pmb)), we find that there are 
several mutations that add [extra positively-charged 
amino acids](${charged_residues}) to the surface. Marine animals typically have these extra charges on 
the surface of their myoglobin 
to help repel neighboring molecules and prevent aggregation when myoglobin is at 
high concentrations.
`;

const description_p2 = `
# Oxygen Bound to Myoglobin
Basic controls for the audio comments :
[Play](${encodeURIComponent(`!play-audio=${_Audio3}`)})
[Pause](!pause-audio)
[Stop](!stop-audio)

A later structure of myoglobin, PDB entry [1mbo](https://www.rcsb.org/structure/1mbo), 
shows that [oxygen](${q('index 1276+1277')}) binds to 
the [iron](${q('index 1275')}) atom deep inside the protein. 
So how does it get in and out? The 
answer is that the structure in the PDB is only one snapshot of the 
protein, caught when it is in a tightly-closed form. In reality, 
myoglobin (and all other proteins) is constantly in motion, performing 
small flexing and breathing motions (illustrated here by PDB entry [1myf](https://www.rcsb.org/structure/1myf)). So, temporary openings constantly 
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
[carbon-rich amino acids](${q('resn ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO')}) are sheltered 
inside and charged amino acids [positively](${q('resn LYS+ARG+HIS')}) 
and [negatively](${q('resn GLU+ASP')}) are most often found on the surface, 
occasionally forming [salt bridges](${q('chain A and resi 44+47+77+18')}) 
that pair two opposite charges (shown here with circles). 

To explore some of these principles, eplxore freely in the interactive view.

# Topics for Further Discussion
You can use the sequence comparison tool to align the sequences of different 
myoglobins, looking for mutations. For instance, [here is the alignment of whale and pig myoglobin](https://www.rcsb.org/alignment?request-body=eJyljrsOwjAMRf%2FFcxgqxNKN7kXsqKpC4pZIeZTEVamq%2FDtOGRAzm%2BPje3I3eM4YV6g3CBOZ4FMZI9IcfZ%2BQoVfYa0kS6kHahFmACp7wReXQBY1QwyRNXExCEOCQHkEX5qUrjNxBWjN64GSiOCtWI%2F9y2wA9xbU3fA1V21w4ndCiKjWKQKbVfeiZ0R3HUmhfVIKz%2Bvs8HXMWv75r2%2Fzn61gJg7F71y6%2FARZEZL8%3D&response-body=eJzlU11vmzAU%2FS88Q2RsjKFv7iCA7GRJQFRVVUU0uClSQjo%2BtkVV%2FvuuyUdTadMeurcJIXF97j33HN%2FLm1HVzzvj5s3o%2B6o0bgzlUaZoyaynwvUshznPll8UpVVgRUr1hAkrmWEabVd0fQv5X75OZjLMQuNgGlvVFZqq2FTreqvqbrndlQqSXouq%2BVG1CgqvMNW97HTLbmsNp5qiUW2%2F6YD44Q16NP2q6%2BFoCKGm2S8HkfbkdqpFqI1addWuHpq2%2B%2B0R5QA9qfWyVd%2BGA9uE2vI9pORwMD%2FyzSa3n%2BN7NN%2FlLi8eB92NWgOl9vDwF1YdVnWpfho3yDQ2ql53L0f%2BR%2FMTtaCta4q6fd4126I7a7FNdHp%2B%2BwUd0chxCXY9xjA2LTRiNkWMONT2TDSimDi%2Bz1yHQDaAmGCbeTbxXR25rodsG%2FvHiCGGEHE9D2vukUcpdgnBlEGAEfU9RrFPdKbDqOMT2x8yLYpH1MWOQxzP9U3CRi5lBCGKoKlFR77j2Rg5iNkgVw%2Bg326LZq%2BH165257X5Xmx62EFo68I97F%2F1PsK99apeKasqYUpVtzf0QlwyXeeSuZikwUfQ9y5gNrGGRoYefw1jr5fQtSp7tdQb3w73rxH9G2xUeUa1MI3Aq2XDEHUpzOxUoKPV7rtqirXSqQjmgdDjcctO0v%2B4ZJ%2FYsQs5lOYyDaPwbi5zGed3XOQhD3IexdE8SGSykGORxrMwk6EYB4uxiKXIQh5OBE%2FDQAoRR3mWy4zLiCcQiiiOAZZiJvk8jXkmYpHMEnEvw3GShjxJYuiTLuJZFIwjHvB5xCdTwWUoxwsRJJyLexHK6H4eDdP453YjmQYnu9P8LrqyG%2BZHu9HFrhjsgs9ru9OT3ejabvbBbn5ld57LeSo%2B2E1PdqfB5Gx3rO3%2BH5t9%2BAU7Kf5z&encoded=true) used to create the illustration in this column.
PDB entry [2jho](https://www.rcsb.org/structure/2jho) includes myoglobin poisoned by cyanide. Take a look and you'll see that the cyanide blocks the binding site for oxygen.


# References

- 1mbn: J. C. Kendrew, R. E. Dickerson, B. E. Strandberg, R. G. Hart, D. R. Davies, D. C. Phillips & V. C. Shore (1960) Structure of Myoglobin. Nature 185, 422-427.
- J. C. Kendrew (1961) The three-dimensional structure of a protein molecule. Scientific American 205(6), 96-110.
- 1mbo: S. E. Phillips (1980) Structure and refinement of oxymyoglobin at 1.6 A resolution. Journal of Molecular Biology 142, 531-554.
- 1pmb: S. J. Smerdon, T. J. Oldfield, E. J. Dodson, G. G. Dodson, R. E. Hubbard & A. J. Wilkinson (1990) Determination of the crystal structure of recombinant pig myoglobin by molecular replacement and its refinement. Acta Crystallographica B, 46, 370-377.
- 1myf: Osapay K, Theriault Y, Wright PE, Case DA. Solution structure of carbonmonoxy myoglobin determined from nuclear magnetic resonance distance and chemical shift constraints. J Mol Biol. 1994;244(2):183-197. doi:10.1006/jmbi.1994.1718
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
            // no outline here

            builder.canvas({ custom: { molstar_postprocessing: { enable_outline: false } } });

            const _1mbn = structure(builder, '1MBN');

            _1mbn.component({ selector: 'ligand' })
                .representation({ ref: 'ligand', type: 'ball_and_stick' })
                .color({ color: 'orange' });
            // FE and O should be spacefill

            _1mbn.component({ selector: { auth_seq_id: 155, label_atom_id: 'F' } })
                .representation({ type: 'spacefill' })
                .color({ color: 'yellow' });

            _1mbn.component({ selector: { auth_seq_id: 154 } })
                .representation({ type: 'spacefill' })
                .color({ color: 'blue' });

            _1mbn.component({ selector: { auth_seq_id: 154 } })
                .representation({ type: 'spacefill' })
                .color({ color: 'blue' });

            const chA = _1mbn.component({ selector: { label_asym_id: 'A' } });
            chA.representation({ type: 'surface', surface_type: 'gaussian' })
                .color({ color: '#ff0303' })
                .opacity({ ref: 'surfopa', opacity: 0.0 });

            chA.representation({ type: 'line' })
                .color({ custom: { molstar_color_theme_name: 'element-symbol' } })
                .opacity({ ref: 'lineopa', opacity: 0.0 });

            chA.representation({ type: 'cartoon' })
                .color({ custom: { molstar_color_theme_name: 'secondary-structure' } });

            // whale
            _1mbn.component({ selector: { label_asym_id: 'A' } })
                .representation({ type: 'spacefill', custom: { molstar_representation_params: { ignoreLight: true } } })
                .colorFromSource({
                    schema: 'all_atomic',
                    category_name: 'atom_site',
                    field_name: 'type_symbol',
                    palette: {
                        kind: 'categorical',
                        colors: {
                            'C': '#FFFFFF',
                            'N': '#CCE6FF',
                            'O': '#FFCCCC',
                            'S': '#FFE680',
                        }
                    }
                }).opacity({ ref: 'cpkopa1', opacity: 0.0 });

            _1mbn.component({ selector: { auth_seq_id: 155 } })
                .representation({ type: 'spacefill', custom: { molstar_representation_params: { ignoreLight: true } } })
                .color({ custom: GColors2 }).opacity({ ref: 'cpkopa2', opacity: 0.0 });

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
                label_size: 5
            });
            addNextButton(builder, 'whale', [13.5, 0, 7.7]);

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
                }
            );
            anim.interpolate({
                kind: 'scalar',
                target_ref: 'lineopa',
                duration_ms: 2000,
                start_ms: 0,
                property: 'opacity',
                start: 0.0,
                end: 1.0,
            });

            anim.interpolate({
                kind: 'scalar',
                target_ref: 'ligand',
                start_ms: 22000,
                duration_ms: 10000,
                frequency: 6,
                alternate_direction: true,
                property: ['custom', 'molstar_representation_params', 'emissive'],
                end: 1.0,
            });

            anim.interpolate({
                kind: 'scalar',
                target_ref: 'cpkopa1',
                duration_ms: 5000,
                start_ms: 40000,
                property: 'opacity',
                start: 0.0,
                end: 1.0,
            });

            anim.interpolate({
                kind: 'scalar',
                target_ref: 'cpkopa2',
                duration_ms: 5000,
                start_ms: 40000,
                property: 'opacity',
                start: 0.0,
                end: 1.0,
            });
            anim.interpolate({
                kind: 'scalar',
                target_ref: 'next',
                duration_ms: 2000,
                start_ms: 43000,
                property: 'label_opacity',
                start: 0.0,
                end: 1.0,
            });
            return builder;
        },
        camera: {
            position: [13.5, 21.1, 73.1],
            target: [13.5, 21.1, 7.7],
            up: [0, 1, 0],
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

            const _1mbn = structure(builder, '1mbn').transform({ ref: 'whalex', translation: [-30, 0, 0] });

            // whale
            _1mbn.component({ selector: { label_asym_id: 'A' } })
                .representation({ type: 'spacefill', custom: { molstar_representation_params: { ignoreLight: true } } })
                .colorFromSource({
                    schema: 'all_atomic', // or maybe just 'atom'
                    category_name: 'atom_site',
                    field_name: 'type_symbol',
                    palette: {
                        kind: 'categorical',
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
                label_tether_length: 1.0,
            })
                .label({
                    text: 'whale',
                    position: { label_asym_id: 'A', auth_seq_id: 8 },
                    label_size: 10
                });

            _1mbn.primitives({
                ref: 'startres',
                label_opacity: 0,
            })
                .label({
                    text: '★', label_offset: 4,
                    position: { label_asym_id: 'A', auth_seq_id: 12, atom_id: 96 }, label_size: 5
                })
                .label({
                    text: '★', label_offset: 4,
                    position: { label_asym_id: 'A', auth_seq_id: 140, auth_atom_id: 'NZ' }, label_size: 5
                })
                .label({
                    text: '★', label_offset: 4,
                    position: { label_asym_id: 'A', auth_seq_id: 87, auth_atom_id: 'NZ' }, label_size: 5
                });

            // the following doesnt work
            const seld = _1mbn.component({
                selector: [
                    { label_asym_id: 'A', auth_seq_id: 12 },
                    { label_asym_id: 'A', auth_seq_id: 140 },
                    { label_asym_id: 'A', auth_seq_id: 87 }
                ]
            });

            seld.representation({ ref: 'scharged', type: 'surface', surface_type: 'gaussian', custom: { molstar_representation_params: { emissive: 0.0, ignoreLight: true } } })
                .colorFromSource(GColors3);

            // pig
            const _1pmb = structure(builder, '1pmb').transform({ ref: 'pig', matrix: align });

            _1pmb.component({ selector: { label_asym_id: 'A' } })
                .representation({ type: 'spacefill', custom: { molstar_representation_params: { ignoreLight: true } } })
                .colorFromSource(GColors3);

            _1pmb.component({ selector: { label_asym_id: 'C', auth_seq_id: 154 } })
                .representation({ type: 'spacefill', custom: { molstar_representation_params: { ignoreLight: true } } })
                .color({ custom: GColors2 });


            _1pmb.primitives({
                ref: 'labelpig',
                label_opacity: 1,
                label_attachment: 'top-center',
                label_show_tether: true,
                label_tether_length: 1.0,
            })
                .label({
                    text: 'pig',
                    position: { label_asym_id: 'A', auth_seq_id: 8 },
                    label_size: 10
                });

            builder.extendRootCustomState({
                molstar_on_load_markdown_commands: {
                    'play-audio': _Audio2,
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
                kind: 'vec3',
                target_ref: 'whalex',
                duration_ms: 10000,
                start_ms: 16000,
                property: 'translation',
                start: [-30, 0, 0],
                end: [-60, 0, 0],
            });
            anim.interpolate({
                kind: 'scalar',
                target_ref: 'startres',
                duration_ms: 1000,
                start_ms: 20000,
                property: 'label_opacity',
                start: 0.0,
                end: 1.0,
            });
            // pig appear at 18s
            anim.interpolate({
                kind: 'transform_matrix',
                target_ref: 'pig',
                duration_ms: 5000,
                start_ms: 18000,
                property: 'matrix',
                translation_start: [-82.54880970106205, 37.49099778180445, -6.133850309914719],
                translation_end: [-52.54880970106205, 37.49099778180445, -6.133850309914719],
            });
            anim.interpolate({
                kind: 'scalar',
                target_ref: 'labelpig',
                duration_ms: 2000,
                start_ms: 18000,
                property: 'label_opacity',
                start: 0.0,
                end: 1.0,
            });
            addNextButton(builder, 'oxygen', [-18.9, 10, 7.3]);
            anim.interpolate({
                kind: 'scalar',
                target_ref: 'next',
                duration_ms: 2000,
                start_ms: 38000,
                property: 'label_opacity',
                start: 0.0,
                end: 1.0,
            });
            anim.interpolate({
                kind: 'scalar',
                target_ref: 'scharged',
                start_ms: 20000,
                duration_ms: 6000,
                frequency: 6,
                alternate_direction: true,
                property: ['custom', 'molstar_representation_params', 'emissive'],
                start: 0.0,
                end: 1.0,
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
            // NMR 1MYF
            // 1A6N unbound
            // 1A6M bound
            // series 2G0R
            const _1mbo = structure(builder, '1mbo')
                .transform({ matrix: alignmbo });

            const _1myf = builder
                .download({ url: pdbUrl('1myf') })
                .parse({ format: 'bcif' })
                .modelStructure({ ref: '1myf' });

            const red1 = '#d3a4a6';
            const red2 = '#d75354';

            const blue1 = '#02d1d1';
            _1myf.component({ selector: { label_asym_id: 'A' } })
                .transform({ translation: [0, 0, 0] })
                .representation({ type: 'spacefill' })
                .color({ color: red1 })
                .opacity({ ref: 'spo', opacity: 1.0 });

            // OXYY
            // should animate in-out in loop
            _1mbo.component({ selector: { label_asym_id: 'C', auth_seq_id: 155 } })
                .representation({ type: 'spacefill' })
                .color({
                    custom: {
                        molstar_color_theme_name: 'element-symbol',
                        molstar_color_theme_params: {
                            carbonColor: {
                                name: 'uniform',
                                params: { value: decodeColor(red2) }
                            },
                        }
                    }
                });

            _1myf.component({ selector: { label_asym_id: 'A' } })
                .representation({ type: 'backbone' })
                .color({ color: red1 });

            _1mbo.component({ selector: { label_asym_id: 'D', auth_seq_id: 555 } })
                .representation({
                    ref: 'oxy', type: 'spacefill', custom: {
                        molstar_representation_params: {
                            emissive: 0.0
                        }
                    }
                })
                .color({ color: blue1 });

            _1mbo.component({ selector: { label_asym_id: 'D', auth_seq_id: 555 } })
                .transform({ ref: 'oxyy', translation: [0, 0, 0] })
                .representation({ type: 'spacefill' })
                .color({ color: blue1 })
                .opacity({ ref: 'oxop', opacity: 0.0 });

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
                property: 'opacity',
                start: 1.0,
                end: 0.05,
            });
            anim.interpolate({
                kind: 'scalar',
                target_ref: '1myf',
                start_ms: 11000,
                duration_ms: 10000,
                frequency: 4,
                alternate_direction: true,
                property: 'model_index',
                start: 0,
                end: 11,
            });
            anim.interpolate({
                kind: 'scalar',
                target_ref: 'oxy',
                start_ms: 3000,
                duration_ms: 10000,
                frequency: 7,
                alternate_direction: true,
                property: ['custom', 'molstar_representation_params', 'emissive'],
                end: 1.0,
            });
            anim.interpolate({
                kind: 'vec3',
                target_ref: 'oxyy',
                duration_ms: 5000,
                start_ms: 16000,
                property: 'translation',
                frequency: 4,
                alternate_direction: false,
                start: [5, -5, -20],
                end: [0, 0, 0],
                noise_magnitude: 1,
            });
            anim.interpolate({
                kind: 'scalar',
                target_ref: 'oxop',
                duration_ms: 1000,
                start_ms: 15000,
                property: 'opacity',
                start: 0.0,
                end: 1.0,
            });
            addNextButton(builder, 'end', [0, -25, 0.0]);
            anim.interpolate({
                kind: 'scalar',
                target_ref: 'next',
                duration_ms: 2000,
                start_ms: 18000,
                property: 'label_opacity',
                start: 0.0,
                end: 1.0,
            });
            return builder;
        },
        camera: {
            position: [-2.2, 0.7, -78.5],
            target: [-0.1, 0.7, 0.6],
            up: [0, 1, 0],
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
            // resn ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO
            const carb = ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO'].map(amk => ({ label_comp_id: amk }));
            // resn LYS+ARG+HIS+ASP+GLU
            const chargedp = ['LYS', 'ARG', 'HIS'].map(amk => ({ label_comp_id: amk }));
            const chargedn = ['ASP', 'GLU'].map(amk => ({ label_comp_id: amk }));

            // salt bridge
            // ASP44-OD1-356-LYS47-NZ-388
            // LYS77-NZ-613-GLU18-OE1-149
            // use primitve distance_measurement
            // and ellipse or ellipsoid with transparancy
            _1mbn.primitives({ ref: 'dist', label_opacity: 0.0 })
                .distance({
                    start: { label_asym_id: 'A', auth_seq_id: 44, atom_id: 356 },
                    end: { label_asym_id: 'A', auth_seq_id: 47, atom_id: 388 },
                    radius: 0.1, dash_length: 0.1,
                    label_size: 2
                })
                .distance({
                    start: { label_asym_id: 'A', auth_seq_id: 77, atom_id: 613 },
                    end: { label_asym_id: 'A', auth_seq_id: 18, atom_id: 149 },
                    radius: 0.1, dash_length: 0.1,
                    label_size: 2
                });
            // 44 OD1 22.300 33.300 -6.200
            // 47 NZ 23.200 32.000 -8.400
            const r44 = Vec3.create(22.300, 33.300, -6.200);
            const r47 = Vec3.create(23.200, 32.000, -8.400);
            getEllipse(builder, r44, r47, 'salt1');

            // 18 OE1 16.600 22.500 20.500
            // 77 NZ 14.100 23.600 22.200
            const r18 = Vec3.create(16.600, 22.500, 20.500);
            const r77 = Vec3.create(14.100, 23.600, 22.200);
            getEllipse(builder, r18, r77, 'salt2');

            const a = _1mbn.component({ selector: carb });
            a.representation({ type: 'ball_and_stick' })
                .color({ color: '#bec0f2' })
                .opacity({ ref: 'carb', opacity: 1.0 });

            const b = _1mbn.component({ selector: chargedp });
            b.representation({ type: 'ball_and_stick' })
                .color({ custom: ill_color('blue', 3.0) })
                .opacity({ ref: 'chargedp', opacity: 1.0 });

            const c = _1mbn.component({ selector: chargedn });
            c.representation({ type: 'ball_and_stick' })
                .color({ custom: ill_color('red', 3.0) })
                .opacity({ ref: 'chargedn', opacity: 1.0 });

            _1mbn.component({ selector: { label_asym_id: 'A' } })
                .representation({ type: 'backbone' })
                .color({ color: '#919191' });

            _1mbn.component({ selector: 'ligand' })
                .representation({
                    ref: 'ligand', type: 'ball_and_stick',
                    custom: {
                        molstar_representation_params: {
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

            anim.interpolate({
                kind: 'scalar',
                target_ref: 'carb',
                duration_ms: 2000,
                start_ms: 8000,
                frequency: 2,
                alternate_direction: true,
                property: 'opacity',
                start: 0.0,
                end: 1.0,
            });

            anim.interpolate({
                kind: 'scalar',
                target_ref: 'chargedp',
                duration_ms: 1000,
                start_ms: 10000,
                property: 'opacity',
                start: 0.0,
                end: 1.0,
            });

            anim.interpolate({
                kind: 'scalar',
                target_ref: 'chargedn',
                duration_ms: 1000,
                start_ms: 10000,
                property: 'opacity',
                start: 0.0,
                end: 1.0,
            });
            // show salt bridge
            anim.interpolate({
                kind: 'scalar',
                target_ref: 'salt1',
                duration_ms: 1000,
                start_ms: 11000,
                property: 'opacity',
                start: 0.0,
                end: 0.3,
            });
            anim.interpolate({
                kind: 'scalar',
                target_ref: 'salt2',
                duration_ms: 1000,
                start_ms: 11000,
                property: 'opacity',
                start: 0.0,
                end: 0.3,
            });
            anim.interpolate({
                kind: 'scalar',
                target_ref: 'dist',
                duration_ms: 1000,
                start_ms: 11000,
                property: 'label_opacity',
                start: 0.0,
                end: 1.0,
            });

            addNextButton(builder, 'intro', [13.5, -10.0, 7.7]);
            anim.interpolate({
                kind: 'scalar',
                target_ref: 'next',
                duration_ms: 2000,
                start_ms: 20000,
                property: 'label_opacity',
                start: 0.0,
                end: 1.0,
            });

            return builder;
        },
        camera: {
            position: [16.0, 47.2, 67.8],
            target: [13.6, 21.1, 7.6],
            up: [0.1, 0.9, -0.4],
        } satisfies MVSNodeParams<'camera'>,
    },
];

function addNextButton(builder: any, snapshotKey: string, position: [number, number, number]) {
    builder.primitives({
        ref: 'next',
        tooltip: 'Click for next part',
        label_opacity: 0,
        label_background_color: 'grey',
        snapshot_key: snapshotKey
    })
        .label({
            ref: 'next_label',
            position: position,
            text: 'Click me to go next',
            label_color: 'white',
            label_size: 5
        });
}
function structure(builder: Root, id: string): MVSStructure {
    return builder
        .download({ url: pdbUrl(id) })
        .parse({ format: 'bcif' })
        .modelStructure();
}

function getEllipse(builder: Root, pos1: Vec3, pos2: Vec3, ref: string) {
    const center = Vec3.add(Vec3(), pos1, pos2);
    Vec3.scale(center, center, 0.5);
    const major_axis = Vec3.sub(Vec3(), pos2, pos1);
    const z_axis = Vec3.create(0, 0, 1);
    // cross to get minor
    const minor_axis = Vec3.cross(Vec3(), major_axis, z_axis);
    return builder.primitives({ ref: ref, opacity: 0.33 }).ellipsoid({
        center: center as any,
        major_axis: major_axis as any,
        minor_axis: minor_axis as any,
        radius: [5.0, 3.0, 3.0],
        color: '#cccccc',
    });
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
            title: 'RCSB Molecule of the Month 1',
            version: '1.0',
            timestamp: new Date().toISOString(),
        }
    };
}