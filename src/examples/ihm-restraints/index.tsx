/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MolViewSpec } from '../../extensions/mvs/behavior';
import { loadMVS } from '../../extensions/mvs/load';
import { createMVSBuilder } from '../../extensions/mvs/tree/mvs/mvs-builder';
import { parseCifText } from '../../mol-io/reader/cif/text/parser';
import { createPluginUI } from '../../mol-plugin-ui';
import { renderReact18 } from '../../mol-plugin-ui/react18';
import { DefaultPluginUISpec } from '../../mol-plugin-ui/spec';
import { PluginConfig } from '../../mol-plugin/config';
import { PluginSpec } from '../../mol-plugin/spec';
import './index.html';
require('../../mol-plugin-ui/skin/light.scss');

async function createViewer(root: HTMLElement) {
    const spec = DefaultPluginUISpec();
    const plugin = await createPluginUI({
        target: root,
        render: renderReact18,
        spec: {
            ...spec,
            layout: {
                initial: {
                    isExpanded: false,
                    showControls: false
                }
            },
            components: {
                remoteState: 'none',
            },
            behaviors: [
                ...spec.behaviors,
                PluginSpec.Behavior(MolViewSpec)
            ],
            config: [
                [PluginConfig.Viewport.ShowAnimation, false],
            ]
        }
    });

    return plugin;
}

async function parseRestraints(url: string) {
    const req = await fetch(url);
    const data = await req.text();

    const parsed = await parseCifText(data).run();

    if (parsed.isError) {
        console.error(parsed);
        return [];
    }

    const dataBlocks = parsed.result.blocks;

    const cat = dataBlocks[0].categories['ihm_cross_link_restraint'];
    const entity_id_1 = cat.getField('entity_id_1');
    const asym_id_1 = cat.getField('asym_id_1');
    const seq_id_1 = cat.getField('seq_id_1');
    const entity_id_2 = cat.getField('entity_id_2');
    const asym_id_2 = cat.getField('asym_id_2');
    const seq_id_2 = cat.getField('seq_id_2');

    const restraints = [];
    for (let i = 0; i < cat.rowCount; i++) {
        restraints.push({
            entity_id_1: entity_id_1?.str(i),
            asym_id_1: asym_id_1?.str(i),
            seq_id_1: seq_id_1?.int(i),
            entity_id_2: entity_id_2?.str(i),
            asym_id_2: asym_id_2?.str(i),
            seq_id_2: seq_id_2?.int(i)
        });
    }

    return restraints;
}

export async function loadIHMRestraints(root: HTMLElement, url?: string) {
    url ??= 'https://pdb-ihm.org/cif/8zz1.cif';

    const plugin = await createViewer(root);
    const restraints = await parseRestraints(url);

    const builder = createMVSBuilder();

    const structure = builder
        .download({ url })
        .parse({ format: 'mmcif' })
        .modelStructure();

    structure
        .component({ selector: 'coarse' })
        .representation({ type: 'spacefill' })
        .color({ custom: { molstar_use_default_coloring: true } });

    structure
        .component({ selector: 'polymer' })
        .representation({ type: 'cartoon' })
        .color({ custom: { molstar_use_default_coloring: true } });

    const primitives = structure.primitives();
    for (const {
        entity_id_1: e1, asym_id_1: a1, seq_id_1: s1,
        entity_id_2: e2, asym_id_2: a2, seq_id_2: s2
    } of restraints) {
        primitives.tube({
            start: { label_entity_id: e1, label_asym_id: a1, label_seq_id: s1 },
            end: { label_entity_id: e2, label_asym_id: a2, label_seq_id: s2 },
            color: 'red',
            radius: 1,
            dash_length: 1,
        });
    }

    const data = builder.getState();

    await loadMVS(plugin, data, { sanityChecks: true, replaceExisting: true });
}

(window as any).loadIHMRestraints = loadIHMRestraints;