/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Column } from '../../../mol-data/db';
import { MmcifFormat } from '../../../mol-model-formats/structure/mmcif';
import { Model } from '../../../mol-model/structure';
import { PluginStateObject } from '../../../mol-plugin-state/objects';
import { StructureComponent } from '../../../mol-plugin-state/transforms/model';
import { StructureRepresentation3D } from '../../../mol-plugin-state/transforms/representation';
import { setSubtreeVisibility } from '../../../mol-plugin/behavior/static/state';
import { PluginContext } from '../../../mol-plugin/context';
import { MolScriptBuilder } from '../../../mol-script/language/builder';
import { ColorNames } from '../../../mol-util/color/names';


/** Amount by which to expand the camera radius when zooming to atoms involved in struct_conn (angstroms) */
const EXTRA_RADIUS = 4;

/** Tags for state tree nodes managed by this extension  */
const TAGS = {
    RESIDUE_SEL: 'structconn-focus-residue-sel',
    ATOM_SEL: 'structconn-focus-atom-sel',
    RESIDUE_REPR: 'structconn-focus-residue-repr',
    RESIDUE_NCI_REPR: 'structconn-focus-residue-nci-repr',
    ATOM_REPR: 'structconn-focus-atom-repr',
} as const;

type VisualParams = ReturnType<typeof StructureRepresentation3D.createDefaultParams>

/** Parameters for 3D representation of atoms involved in struct_conn (pink bubbles) */
const ATOMS_VISUAL_PARAMS: VisualParams = {
    type: { name: 'ball-and-stick', params: { sizeFactor: 0.25, sizeAspectRatio: 0.73, adjustCylinderLength: true, xrayShaded: true, aromaticBonds: false, multipleBonds: 'off', dashCount: 1, dashCap: false } },
    colorTheme: { name: 'uniform', params: { value: ColorNames.magenta } },
    sizeTheme: { name: 'physical', params: {} },
} as const;

/** Parameters for 3D representation of residues involved in struct_conn (normal ball-and-stick) */
const RESIDUES_VISUAL_PARAMS: VisualParams = {
    type: { name: 'ball-and-stick', params: { sizeFactor: 0.16 } },
    colorTheme: { name: 'element-symbol', params: {} },
    sizeTheme: { name: 'physical', params: {} },
} as const;


/** All public functions provided by the StructConn extension  */
export const wwPDBStructConnExtensionFunctions = {
    /** Return an object with all struct_conn records for a loaded structure.
     * Applies to the first structure belonging to `entry` (e.g. '1tqn'),
     * or to the first loaded structure overall if `entry` is `undefined`.
     */
    getStructConns(plugin: PluginContext, entry: string | undefined): { [id: string]: StructConnRecord } {
        const structNode = selectStructureNode(plugin, entry);
        const structure = structNode?.obj?.data;
        if (structure) return extractStructConns(structure.model);
        else return {};
    },

    /** Create visuals for residues and atoms involved in a struct_conn with ID `structConnId`
     * and zoom on them. If `keepExisting` is false (default), remove any such visuals created by previous calls to this function.
     * Also hide all carbohydrate SNFG visuals within the structure (as they would occlude our residues of interest).
     * Return a promise that resolves to the number of involved atoms which were successfully selected (2, 1, or 0).
     */
    async inspectStructConn(plugin: PluginContext, entry: string | undefined, structConnId: string, keepExisting: boolean = false): Promise<number> {
        const structNode = selectStructureNode(plugin, entry);
        const structure = structNode?.obj?.data;
        if (!structure) {
            console.error('Structure not found');
            return 0;
        }

        const conns: { [id: string]: StructConnRecord } = structure.model._staticPropertyData['wwpdb-struct-conn-extension-data'] ??= extractStructConns(structure.model);
        const conn = conns[structConnId];
        if (!conn) {
            console.error(`The structure does not contain struct_conn "${structConnId}"`);
            return 0;
        }

        if (!keepExisting) {
            await removeAllStructConnInspections(plugin, structNode);
        }
        const nSelectedAtoms = await addStructConnInspection(plugin, structNode, conn);
        hideSnfgNodes(plugin, structNode);
        return nSelectedAtoms;
    },

    /** Remove anything created by `inspectStructConn` within the structure and
     * make visible any carbohydrate SNFG visuals that have been hidden by `inspectStructConn`.
     */
    async clearStructConnInspections(plugin: PluginContext, entry: string | undefined) {
        const structNode = selectStructureNode(plugin, entry);
        if (!structNode) return;
        await removeAllStructConnInspections(plugin, structNode);
        unhideSnfgNodes(plugin, structNode);
    },
};


type StructNode = Exclude<ReturnType<typeof selectStructureNode>, undefined>

/** Return the first structure node belonging to `entry` (e.g. '1tqn'),
 * or to the first loaded structure node overall if `entry` is `undefined`.
 * Includes only "root" structures, not structure components. */
function selectStructureNode(plugin: PluginContext, entry: string | undefined) {
    const structNodes = plugin.state.data
        .selectQ(q => q.rootsOfType(PluginStateObject.Molecule.Structure));
    if (entry) {
        const result = structNodes.find(node => node.obj && node.obj.data.model.entry.toLowerCase() === entry.toLowerCase());
        if (!result) {
            console.warn(`Structure with entry ID "${entry}" was not found. Available structures: ${structNodes.map(node => node.obj?.data.model.entry)}`);
        }
        return result;
    } else {
        if (structNodes.length > 1) {
            console.warn(`Structure entry ID was not specified, but there is more than one loaded structure (${structNodes.map(node => node.obj?.data.model.entry)}). Taking the first structure.`);
        }
        if (structNodes.length === 0) {
            console.warn(`There are no loaded structures.`);
        }
        return structNodes[0];
    }
}


/** Represents one partner (i.e. atom) of a struct_conn */
interface StructConnPartner {
    asymId: string,
    seqId: number | undefined,
    authSeqId: number | undefined,
    insCode: string,
    compId: string,
    atomId: string,
    /** Alternative location (use empty string if not given) */
    altId: string,
}

/** Represents a struct_conn (interaction between two partners) */
export interface StructConnRecord {
    id: string,
    distance: number,
    partner1: StructConnPartner,
    partner2: StructConnPartner,
}


/** Return an object with all struct_conn records read from mmCIF.
 * Return {} if the model comes from another format than mmCIF.
 */
function extractStructConns(model: Model): { [id: string]: StructConnRecord } {
    if (!MmcifFormat.is(model.sourceData)) {
        console.error('Cannot get struct_conn because source data are not mmCIF.');
        return {};
    }
    const mmcifData = model.sourceData.data;
    const {
        id,
        ptnr1_label_asym_id: asym1,
        ptnr1_label_seq_id: seq1,
        ptnr1_auth_seq_id: authSeq1,
        pdbx_ptnr1_PDB_ins_code: authInsCode1,
        ptnr1_label_comp_id: comp1,
        ptnr1_label_atom_id: atom1,
        pdbx_ptnr1_label_alt_id: alt1,
        ptnr2_label_asym_id: asym2,
        ptnr2_label_seq_id: seq2,
        ptnr2_auth_seq_id: authSeq2,
        pdbx_ptnr2_PDB_ins_code: authInsCode2,
        ptnr2_label_comp_id: comp2,
        ptnr2_label_atom_id: atom2,
        pdbx_ptnr2_label_alt_id: alt2,
        pdbx_dist_value: distance } = mmcifData.db.struct_conn;
    const n = id.rowCount;
    const result: { [id: string]: StructConnRecord } = {};
    for (let i = 0; i < n; i++) {
        const conn: StructConnRecord = {
            id: id.value(i),
            distance: distance.value(i),
            partner1: {
                asymId: asym1.value(i),
                seqId: seq1.valueKind(i) === Column.ValueKinds.Present ? seq1.value(i) : undefined,
                authSeqId: authSeq1.valueKind(i) === Column.ValueKinds.Present ? authSeq1.value(i) : undefined,
                insCode: authInsCode1.value(i),
                compId: comp1.value(i),
                atomId: atom1.value(i),
                altId: alt1.value(i),
            },
            partner2: {
                asymId: asym2.value(i),
                seqId: seq2.valueKind(i) === Column.ValueKinds.Present ? seq2.value(i) : undefined,
                authSeqId: authSeq2.valueKind(i) === Column.ValueKinds.Present ? authSeq2.value(i) : undefined,
                insCode: authInsCode2.value(i),
                compId: comp2.value(i),
                atomId: atom2.value(i),
                altId: alt2.value(i),
            },
        };
        result[conn.id] = conn;
    }
    return result;
}

/** Return MolScript expression for atoms or residues involved in a struct_conn */
function structConnExpression(conn: StructConnRecord, by: 'atoms' | 'residues') {
    const { core, struct } = MolScriptBuilder;
    const partnerExpressions = [];
    for (const partner of [conn.partner1, conn.partner2]) {
        const propTests: Parameters<typeof struct.generator.atomGroups>[0] = {
            'chain-test': core.rel.eq([struct.atomProperty.macromolecular.label_asym_id(), partner.asymId]),
            'group-by': struct.atomProperty.core.operatorName(),
        };
        if (partner.seqId !== undefined) {
            propTests['residue-test'] = core.rel.eq([struct.atomProperty.macromolecular.label_seq_id(), partner.seqId]);
        } else if (partner.authSeqId !== undefined) { // for the case of water and carbohydrates (see 5elb, covale3 vs covale5)
            propTests['residue-test'] = core.logic.and([
                core.rel.eq([struct.atomProperty.macromolecular.auth_seq_id(), partner.authSeqId]),
                core.rel.eq([struct.atomProperty.macromolecular.pdbx_PDB_ins_code(), partner.insCode]),
            ]);
        }
        if (by === 'residues' && partner.altId !== '') {
            propTests['atom-test'] = core.rel.eq([struct.atomProperty.macromolecular.label_alt_id(), partner.altId]);
        }
        if (by === 'atoms') {
            propTests['atom-test'] = core.logic.and([
                core.rel.eq([struct.atomProperty.macromolecular.label_atom_id(), partner.atomId]),
                core.rel.eq([struct.atomProperty.macromolecular.label_alt_id(), partner.altId]),
            ]);
        }
        partnerExpressions.push(struct.filter.first([struct.generator.atomGroups(propTests)]));
    }
    return struct.combinator.merge(partnerExpressions.map(e => struct.modifier.union([e])));
}

/** Create visuals for residues and atoms involved in a struct_conn and zoom on them.
 * Return a promise that resolves to the number of involved atoms which were successfully selected (2, 1, or 0).
 */
async function addStructConnInspection(plugin: PluginContext, structNode: StructNode, conn: StructConnRecord): Promise<number> {
    const expressionByResidues = structConnExpression(conn, 'residues');
    const expressionByAtoms = structConnExpression(conn, 'atoms');

    const update = plugin.build();
    update.to(structNode).apply(
        StructureComponent,
        { label: `${conn.id} (residues)`, type: { name: 'expression', params: expressionByResidues } },
        { tags: [TAGS.RESIDUE_SEL] }
    ).apply(
        StructureRepresentation3D,
        RESIDUES_VISUAL_PARAMS,
        { tags: [TAGS.RESIDUE_REPR] }
    );

    const atomsSelection = update.to(structNode).apply(
        StructureComponent,
        { label: `${conn.id} (atoms)`, type: { name: 'expression', params: expressionByAtoms } },
        { tags: [TAGS.ATOM_SEL] }
    );
    const atomsVisual = update.to(atomsSelection.ref).apply(
        StructureRepresentation3D,
        ATOMS_VISUAL_PARAMS,
        { tags: [TAGS.ATOM_REPR] }
    );
    await update.commit();

    plugin.managers.camera.focusRenderObjects(atomsVisual.selector.data?.repr.renderObjects, { extraRadius: EXTRA_RADIUS });
    const nSelectedAtoms = atomsSelection.selector.obj?.data?.elementCount ?? 0;
    return nSelectedAtoms;
}

/** Remove anything created by `addStructConnInspection` */
async function removeAllStructConnInspections(plugin: PluginContext, structNode: StructNode) {
    const selNodes = [
        ...plugin.state.data.selectQ(q => q.byRef(structNode.transform.ref).subtree().withTag(TAGS.RESIDUE_SEL)),
        ...plugin.state.data.selectQ(q => q.byRef(structNode.transform.ref).subtree().withTag(TAGS.ATOM_SEL)),
    ];
    const update = plugin.build();
    for (const node of selNodes) {
        update.delete(node);
    }
    await update.commit();
}

/** Hide all carbohydrate SNFG visuals */
function hideSnfgNodes(plugin: PluginContext, structNode: StructNode) {
    const snfgNodes = plugin.state.data.selectQ(q => q.byRef(structNode.transform.ref).subtree().withTag('branched-snfg-3d'));
    for (const node of snfgNodes) {
        setSubtreeVisibility(plugin.state.data, node.transform.ref, true); // true means hidden
    }
}

/** Make visible all carbohydrate SNFG visuals that have been hidden by `hideSnfgNodes` */
function unhideSnfgNodes(plugin: PluginContext, structNode: StructNode) {
    const snfgNodes = plugin.state.data.selectQ(q => q.byRef(structNode.transform.ref).subtree().withTag('branched-snfg-3d'));
    for (const node of snfgNodes) {
        try {
            setSubtreeVisibility(plugin.state.data, node.transform.ref, false); // false means visible
        } catch {
            // this is OK, the node has been removed
        }
    }
}
