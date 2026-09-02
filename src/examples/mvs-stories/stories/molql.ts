/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MVSData_States } from '../../../extensions/mvs/mvs-data';
import { createMVSBuilder } from '../../../extensions/mvs/tree/mvs/mvs-builder';
import { MolScriptBuilder as MS } from '../../../mol-script/language/builder';
import { parse } from '../../../mol-script/transpile';

const imatinib = MS.struct.generator.atomGroups({
    'chain-test': MS.core.rel.eq([
        MS.struct.atomProperty.macromolecular.label_asym_id(),
        'G',
    ]),
});

const gatekeeper = MS.struct.generator.atomGroups({
    'chain-test': MS.core.rel.eq([
        MS.struct.atomProperty.macromolecular.auth_asym_id(),
        'A',
    ]),
    'residue-test': MS.core.rel.eq([
        MS.struct.atomProperty.macromolecular.auth_seq_id(),
        315,
    ]),
});

const imatinibN13 = MS.struct.generator.atomGroups({
    'chain-test': MS.core.rel.eq([
        MS.struct.atomProperty.macromolecular.label_asym_id(),
        'G',
    ]),
    'atom-test': MS.core.rel.eq([
        MS.struct.atomProperty.macromolecular.label_atom_id(),
        'N13',
    ]),
});

// Transpilers resolve their input syntax to the same JSON MolQL tree as the builder above.
const bindingPocket = parse('pymol', 'byres polymer within 5 of resn STI');
const thr315OG1 = parse('pymol', 'chain A and resi 315 and name OG1');

export function buildStory(): MVSData_States {
    return {
        kind: 'multiple',
        snapshots: [componentSnapshot(), colorSnapshot(), interactionSnapshot(), primitiveSnapshot()],
        metadata: {
            title: 'Experimental MolQL Selectors',
            version: '1.0',
            timestamp: new Date().toISOString(),
        },
    };
}

function componentSnapshot() {
    const builder = createMVSBuilder();
    const structure = loadStructure(builder);

    structure.component({ selector: 'polymer' })
        .representation({ type: 'cartoon' })
        .color({ color: '#8AA6C1' });

    structure.component({ selector: { molql: imatinib } })
        .representation({ type: 'ball_and_stick' })
        .color({ color: '#F08A4B' });

    structure.component({ selector: { molql: bindingPocket } })
        .representation({ type: 'ball_and_stick' })
        .color({ color: '#B8497A' });

    return builder.getSnapshot({
        title: 'MolQL component selectors',
        key: 'components',
        description: `# MolQL component selectors

Every non-static component selector is wrapped as \`{ molql: ... }\`. The orange imatinib selection is built with \`MolScriptBuilder\`; the pink binding pocket is compiled from the PyMOL expression \`byres polymer within 5 of resn STI\`.

Both are stored in MolViewSpec as JSON MolQL expression trees, not executable source text.`,
        description_format: 'markdown',
    });
}

function colorSnapshot() {
    const builder = createMVSBuilder();
    const structure = loadStructure(builder);
    const representation = structure.component({ selector: 'polymer' })
        .representation({ type: 'cartoon' });

    representation.color({ color: '#8AA6C1' });
    representation.color({ selector: { molql: bindingPocket }, color: '#B8497A' });
    representation.color({ selector: { molql: gatekeeper }, color: '#F08A4B' });

    return builder.getSnapshot({
        title: 'MolQL color selectors',
        key: 'colors',
        description: `# MolQL color selectors

MolQL is also accepted by the \`selector\` on a \`color\` node. Here one broad cartoon representation is selectively recolored: the PyMOL-transpiled pocket is pink and the builder-generated gatekeeper residue Thr315 is orange.

This demonstrates multilayer coloring without creating separate visual components for each selected region.`,
        description_format: 'markdown',
    });
}

function interactionSnapshot() {
    const builder = createMVSBuilder();
    const structure = loadStructure(builder);

    structure.component({ selector: 'polymer' })
        .representation({ type: 'cartoon' })
        .color({ color: '#8AA6C1' });

    structure.component({ selector: { molql: imatinib } })
        .representation({ type: 'ball_and_stick' })
        .color({ color: '#F08A4B' });

    const pocket = structure.component({ selector: { molql: bindingPocket } });
    pocket.representation({ type: 'ball_and_stick' }).color({ color: '#B8497A' });
    pocket.label({ text: 'PyMOL binding-pocket selection' });
    pocket.tooltip({ text: 'Polymer residues within 5 Å of imatinib (STI).' });

    return builder.getSnapshot({
        title: 'MolQL component interactions',
        key: 'interactions',
        description: `# Reusing a MolQL component

Once a wrapped MolQL selector creates a component, its child labels and tooltips operate on that computed component as usual. Hover the pink pocket to see its tooltip.

This keeps the query in the component node while labels, representations, and other component children reuse its result.`,
        description_format: 'markdown',
    });
}

function primitiveSnapshot() {
    const builder = createMVSBuilder();
    const structure = loadStructure(builder);

    structure.component({ selector: 'polymer' })
        .representation({ type: 'cartoon' })
        .color({ color: '#8AA6C1' });
    structure.component({ selector: { molql: imatinib } })
        .representation({ type: 'ball_and_stick' })
        .color({ color: '#F08A4B' });
    structure.component({ selector: { molql: gatekeeper } })
        .representation({ type: 'ball_and_stick' })
        .color({ color: '#B8497A' });

    const primitives = structure.primitives();
    primitives.distance({
        start: { molql: imatinibN13 },
        end: { molql: thr315OG1 },
        color: '#F08A4B',
        radius: 0.12,
        dash_length: 0.2,
        label_template: 'Imatinib N13–Thr315 OG1: {{distance}}',
        label_color: '#F08A4B',
    });

    return builder.getSnapshot({
        title: 'MolQL primitive positions',
        key: 'primitives',
        description: `# MolQL primitive positions

Primitive positions accept wrapped MolQL too. This dashed distance measurement connects imatinib atom N13, selected with \`MolScriptBuilder\`, to Thr315 atom OG1, selected from the PyMOL expression \`chain A and resi 315 and name OG1\`.

Because both queries select one atom, their boundary-sphere centers are the exact atom coordinates. This makes the primitive a real measurement of a known imatinib–Thr315 hydrogen-bond contact.`,
        description_format: 'markdown',
    });
}

function loadStructure(builder: ReturnType<typeof createMVSBuilder>) {
    return builder
        .download({ url: 'https://www.ebi.ac.uk/pdbe/entry-files/download/1iep.bcif' })
        .parse({ format: 'bcif' })
        .assemblyStructure();
}
