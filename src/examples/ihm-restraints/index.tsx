/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MolViewSpec } from '../../extensions/mvs/behavior';
import { loadMVS } from '../../extensions/mvs/load';
import { MVSData_States, Snapshot } from '../../extensions/mvs/mvs-data';
import { createMVSBuilder } from '../../extensions/mvs/tree/mvs/mvs-builder';
import { parseCifText } from '../../mol-io/reader/cif/text/parser';
import { Vec3 } from '../../mol-math/linear-algebra';
import { trajectoryFromMmCIF } from '../../mol-model-formats/structure/mmcif';
import { Model } from '../../mol-model/structure';
import { CoarseElementKey, CoarseElementReference } from '../../mol-model/structure/model/properties/coarse';
import { createPluginUI } from '../../mol-plugin-ui';
import { renderReact18 } from '../../mol-plugin-ui/react18';
import { DefaultPluginUISpec } from '../../mol-plugin-ui/spec';
import { PluginConfig } from '../../mol-plugin/config';
import { PluginContext } from '../../mol-plugin/context';
import { PluginSpec } from '../../mol-plugin/spec';
import { Task } from '../../mol-task';
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
                    isExpanded: true,
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
                [PluginConfig.Viewport.ShowTrajectoryControls, false],
            ]
        }
    });

    return plugin;
}

interface IHMRestraintInfo {
    e1: CoarseElementKey & { label_comp_id: string },
    e2: CoarseElementKey & { label_comp_id: string },
    a: Vec3,
    b: Vec3,
    restraintType: 'harmonic' | 'upper bound' | 'lower bound',
    threshold: number,
    satisfied: boolean,
    distance: number,
}


function getCoarseElementPosition(e: CoarseElementReference, model: Model, position: Vec3) {
    if (!e.kind) Vec3.set(position, 0, 0, 0);
    const { x, y, z } = model.coarseConformation[e.kind!];
    const idx = e.index;
    Vec3.set(position, x[idx], y[idx], z[idx]);
}

const _elementRef = CoarseElementReference();
function resolvePosition(model: Model, key: CoarseElementKey, position: Vec3) {
    if (model.coarseHierarchy.index.findElement(key, _elementRef)) {
        getCoarseElementPosition(_elementRef, model, position);
        return true;
    }

    const rI = model.atomicHierarchy.index.findResidueLabel(key);
    if (rI < 0) return false;

    const atomStart = model.atomicHierarchy.residueAtomSegments.offsets[rI];
    const atomEnd = model.atomicHierarchy.residueAtomSegments.offsets[rI + 1];
    const atomId = model.atomicHierarchy.atoms.label_atom_id;
    let aI = atomStart;
    // Find CA otherwise use the first atom.
    // Possible future improvement: use the atom closest to the center of mass of the residue.
    for (; aI < atomEnd; aI++) {
        if (atomId.value(aI) === 'CA') break;
    }
    if (aI === atomEnd) aI = atomStart;

    const { x, y, z } = model.atomicConformation;
    Vec3.set(position, x[aI], y[aI], z[aI]);

    return true;
}

const HarmonicRestraintTolerance = 0.1;

async function parseRestraints(plugin: PluginContext, url: string) {
    const req = await fetch(url);
    const data = await req.text();

    const parsed = await plugin.runTask(parseCifText(data), { useOverlay: true });

    if (parsed.isError) {
        console.error(parsed);
        return [];
    }

    const trajectory = await plugin.runTask(trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result));

    const dataBlocks = parsed.result.blocks;

    const cat = dataBlocks[0].categories['ihm_cross_link_restraint'];
    const entity_id_1 = cat.getField('entity_id_1')!;
    const asym_id_1 = cat.getField('asym_id_1')!;
    const seq_id_1 = cat.getField('seq_id_1')!;
    const comp_id_1 = cat.getField('comp_id_1')!;
    const entity_id_2 = cat.getField('entity_id_2')!;
    const asym_id_2 = cat.getField('asym_id_2')!;
    const seq_id_2 = cat.getField('seq_id_2')!;
    const comp_id_2 = cat.getField('comp_id_2')!;
    const restraint_type = cat.getField('restraint_type')!;
    const threshold = cat.getField('distance_threshold')!;

    const e1key = CoarseElementKey();
    const e2key = CoarseElementKey();

    const a = Vec3.zero();
    const b = Vec3.zero();

    const modelRestraints: IHMRestraintInfo[][] = [];

    for (let modelIndex = 0; modelIndex < trajectory.frameCount; modelIndex++) {
        const _model = trajectory.getFrameAtIndex(modelIndex);
        const model = Task.is(_model) ? await plugin.runTask(_model) : _model;

        const restraints: IHMRestraintInfo[] = [];
        modelRestraints.push(restraints);

        for (let i = 0; i < cat.rowCount; i++) {
            e1key.label_entity_id = entity_id_1.str(i);
            e1key.label_asym_id = asym_id_1.str(i);
            e1key.label_seq_id = seq_id_1.int(i);
            e2key.label_entity_id = entity_id_2.str(i);
            e2key.label_asym_id = asym_id_2.str(i);
            e2key.label_seq_id = seq_id_2.int(i);

            if (!resolvePosition(model, e1key, a) || !resolvePosition(model, e2key, b)) {
                continue;
            }

            const restraintType: 'harmonic' | 'upper bound' | 'lower bound' = restraint_type.str(i)?.toLowerCase() as any;
            const thresholdValue = threshold.float(i);
            const distance = Vec3.distance(a, b);

            let satisfied = true;
            if (restraintType === 'harmonic') {
                const thresholdValue = threshold.float(i);
                satisfied = distance >= (1 - HarmonicRestraintTolerance) * thresholdValue && distance <= (1 + HarmonicRestraintTolerance) * thresholdValue;
            } else if (restraintType === 'upper bound') {
                satisfied = distance <= thresholdValue;
            } else if (restraintType === 'lower bound') {
                satisfied = distance >= thresholdValue;
            }

            restraints.push({
                e1: { ...e1key, label_comp_id: comp_id_1.str(i) },
                e2: { ...e2key, label_comp_id: comp_id_2.str(i) },
                a: Vec3.clone(a),
                b: Vec3.clone(b),
                restraintType,
                threshold: thresholdValue,
                satisfied,
                distance,
            });
        }
    }

    return modelRestraints;
}

function baseStructure(url: string, modelIndex: number) {
    const builder = createMVSBuilder();

    const structure = builder
        .download({ url })
        .parse({ format: 'mmcif' })
        .modelStructure({ model_index: modelIndex });

    structure
        .component({ selector: 'coarse' })
        .representation({ type: 'spacefill' })
        .color({ custom: { molstar_use_default_coloring: true } })
        .opacity({ opacity: 0.51 });

    structure
        .component({ selector: 'polymer' })
        .representation({ type: 'cartoon' })
        .color({ custom: { molstar_use_default_coloring: true } })
        .opacity({ opacity: 0.51 });

    return [builder, structure] as const;
}

function drawConstraints([, structure]: ReturnType<typeof baseStructure>, restraints: IHMRestraintInfo[], options: {
    filter: (r: IHMRestraintInfo) => boolean,
    color: (r: IHMRestraintInfo) => any,
    radius?: (r: IHMRestraintInfo) => number,
    tooltip: (r: IHMRestraintInfo) => string | undefined,
}) {
    const primitives = structure.primitives();
    for (const r of restraints) {
        if (!options.filter(r)) continue;

        const radius = options.radius?.(r) ?? 1;

        primitives.tube({
            start: r.a as any,
            end: r.b as any,
            color: options.color(r) || 'white',
            tooltip: options.tooltip(r),
            radius: radius,
            dash_length: radius,
        });
    }
}

function restraintTooltip(r: IHMRestraintInfo) {
    return `
- Distance: ${r.distance.toFixed(2)} Å
- Threshold: ${r.threshold.toFixed(2)} Å
- Constraint: ${r.restraintType}
- Satisfied: ${r.satisfied ? 'Yes' : 'No'}
`;
}

export async function loadIHMRestraints(root: HTMLElement, url?: string) {
    url ??= 'https://pdb-ihm.org/cif/8zz1.cif';

    const plugin = await createViewer(root);
    const modelRestraints = await parseRestraints(plugin, url);

    const modelIndex = 0;
    const restraints = modelRestraints[modelIndex];

    const nVialoted = restraints.filter(r => !r.satisfied).length;
    const nSatisfied = restraints.length - nVialoted;

    const snapshots: Snapshot[] = [];

    let mvs = baseStructure(url, modelIndex);
    drawConstraints(mvs, restraints, {
        filter: r => true,
        color: r => r.e1.label_entity_id === r.e2.label_entity_id && r.e1.label_asym_id === r.e2.label_asym_id ? 'yellow' : 'blue',
        radius: r => 1,
        tooltip: restraintTooltip,
    });
    snapshots.push(mvs[0].getSnapshot({
        title: 'All Restraints',
        linger_duration_ms: 5000,
        description: `
### All Restraints

- Yellow: Intra-chain restraints
- Blue: Inter-chain restraints
`,
    }));

    mvs = baseStructure(url, modelIndex);
    drawConstraints(mvs, restraints, {
        filter: r => true,
        color: r => r.satisfied ? 'green' : 'red',
        radius: r => 1,
        tooltip: restraintTooltip,
    });
    snapshots.push(mvs[0].getSnapshot({
        title: 'Restraint Validation',
        linger_duration_ms: 5000,
        description: `
### Restraint Validation

- Red: ${nVialoted} Violated restraints
- Green: ${nSatisfied} Satisfied restraints
`,
    }));

    mvs = baseStructure(url, modelIndex);
    drawConstraints(mvs, restraints, {
        filter: r => !r.satisfied,
        color: r => r.satisfied ? 'green' : 'red',
        radius: r => 1,
        tooltip: restraintTooltip,
    });
    snapshots.push(mvs[0].getSnapshot({
        title: 'Violated Restraints',
        linger_duration_ms: 5000,
        description: `
### Violated Restraints

${nVialoted} restraints are violated.
`,
    }));

    mvs = baseStructure(url, modelIndex);
    drawConstraints(mvs, restraints, {
        filter: r => r.satisfied,
        color: r => r.satisfied ? 'green' : 'red',
        radius: r => 1,
        tooltip: restraintTooltip,
    });
    snapshots.push(mvs[0].getSnapshot({
        title: 'Satisfied Restraints',
        linger_duration_ms: 5000,
        description: `
### Violated Restraints

${nSatisfied} restraints are violated.
`,
    }));

    const data: MVSData_States = {
        kind: 'multiple',
        snapshots,
        metadata: {
            title: 'I/HM Restraints',
            version: '1.0',
            timestamp: new Date().toISOString(),
        }
    };

    await loadMVS(plugin, data, { sanityChecks: true, replaceExisting: true, keepSnapshotCamera: true });
}

(window as any).loadIHMRestraints = loadIHMRestraints;