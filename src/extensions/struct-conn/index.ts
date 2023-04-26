/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Column } from '../../mol-data/db';
import { MmcifFormat } from '../../mol-model-formats/structure/mmcif';
import { Structure } from '../../mol-model/structure';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { StructureComponent } from '../../mol-plugin-state/transforms/model';
import { StructureRepresentation3D } from '../../mol-plugin-state/transforms/representation';
import { PluginBehavior } from '../../mol-plugin/behavior';
import { setSubtreeVisibility } from '../../mol-plugin/behavior/static/state';
import { PluginCommand } from '../../mol-plugin/command';
import { PluginContext } from '../../mol-plugin/context';
import { MolScriptBuilder } from '../../mol-script/language/builder';
import { ColorNames } from '../../mol-util/color/names';


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
export const StructConnExtensionFunctions = {
    /** Return an object with all struct_conn records for a loaded structure.
     * Applies to the first structure belonging to `entry` (e.g. '1tqn'),
     * or to the first loaded structure overall if `entry` is `undefined`.
     */
    getStructConns(plugin: PluginContext, entry: string | undefined): { [id: string]: StructConnRecord } {
        const structNode = selectStructureNode(plugin, entry);
        const structure = structNode?.obj?.data;
        if (structure) return extractStructConns(structure);
        else return {};
    },

    /** Create visuals for residues and atoms involved in a struct_conn with ID `structConnId`
     * and zoom on them. If `keepExisting` is false (default), remove any such visuals created by previous calls to this function.
     * Also hide all carbohydrate SNFG visuals (as they would occlude our residues of interest).
     */
    async inspectStructConn(plugin: PluginContext, entry: string | undefined, structConnId: string, keepExisting: boolean = false) {
        const structNode = selectStructureNode(plugin, entry);
        const structure = structNode?.obj?.data;
        if (!structure) {
            console.error('Structure not found');
            return;
        }

        const structNodeRef = structNode.transform.ref;
        const structConnCache = getExtensionState(plugin).structConnCache;
        if (!structConnCache.has(structNodeRef)) {
            structConnCache.set(structNodeRef, extractStructConns(structure));
        }
        const conns = structConnCache.get(structNodeRef)!;
        const conn = conns[structConnId];
        if (!conn) {
            console.error(`The structure does not contain struct_conn "${structConnId}"`);
            return;
        }

        if (!keepExisting) {
            await removeAllStructConnInspections(plugin);
        }
        await addStructConnInspection(plugin, structNode, conn);
        hideSnfgNodes(plugin);
    },

    /** Remove anything created by `inspectStructConn` and
     * make visible any carbohydrate SNFG visuals that have been hidden by `inspectStructConn`.
     */
    async clearStructConnInspections(plugin: PluginContext) {
        await removeAllStructConnInspections(plugin);
        unhideSnfgNodes(plugin);
    },
};


type StructNode = ReturnType<typeof selectStructureNode>

/** Return the first structure node belonging to `entry` (e.g. '1tqn'),
 * or to the first loaded structure node overall if `entry` is `undefined`.
 * Includes only "root" structures, not structure components. */
function selectStructureNode(plugin: PluginContext, entry: string | undefined) {
    const structNodes = plugin.state.data
        .selectQ(q => q.ofType(PluginStateObject.Molecule.Structure))
        .filter(node => node.obj && !node.obj.data.parent && !node.transform.transformer.definition.isDecorator);
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
 * Return {} if structure comes from another format than mmCIF.
 */
function extractStructConns(structure: Structure): { [id: string]: StructConnRecord } {
    if (!MmcifFormat.is(structure.model.sourceData)) {
        console.error('Cannot get struct_conn because source data are not mmCIF.');
        return {};
    }
    const mmcifData = structure.model.sourceData.data;
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

/** Create visuals for residues and atoms involved in a struct_conn and zoom on them */
async function addStructConnInspection(plugin: PluginContext, structNode: StructNode, conn: StructConnRecord) {
    if (!structNode) return;

    const expressionByResidues = structConnExpression(conn, 'residues');
    const expressionByAtoms = structConnExpression(conn, 'atoms');

    const selectionByResidues = await plugin.build().to(structNode).apply(
        StructureComponent,
        { label: `${conn.id} (residues)`, type: { name: 'expression', params: expressionByResidues } },
        { tags: [TAGS.RESIDUE_SEL] }
    ).commit();
    const selectionByAtom = await plugin.build().to(structNode).apply(
        StructureComponent,
        { label: `${conn.id} (atoms)`, type: { name: 'expression', params: expressionByAtoms } },
        { tags: [TAGS.ATOM_SEL] }
    ).commit();

    if (selectionByResidues.isOk) {
        await plugin.build().to(selectionByResidues).apply(
            StructureRepresentation3D,
            RESIDUES_VISUAL_PARAMS,
            { tags: [TAGS.RESIDUE_REPR] }
        ).commit();
    }
    if (selectionByAtom.isOk) {
        const atomBalls = await plugin.build().to(selectionByAtom).apply(
            StructureRepresentation3D,
            ATOMS_VISUAL_PARAMS,
            { tags: [TAGS.ATOM_REPR] }
        ).commit();
        plugin.managers.camera.focusRenderObjects(atomBalls.data?.repr.renderObjects, { extraRadius: EXTRA_RADIUS });
    }
}

/** Remove anything created by `addStructConnInspection` */
async function removeAllStructConnInspections(plugin: PluginContext) {
    const selNodes = [
        ...plugin.state.data.selectQ(q => q.root.subtree().withTag(TAGS.RESIDUE_SEL)),
        ...plugin.state.data.selectQ(q => q.root.subtree().withTag(TAGS.ATOM_SEL)),
    ];
    for (const node of selNodes) {
        await plugin.build().delete(node).commit();
    }
}

/** Hide all carbohydrate SNFG visuals */
function hideSnfgNodes(plugin: PluginContext) {
    const hiddenNodes = getExtensionState(plugin).hiddenNodes;
    const snfgNodes = plugin.state.data.selectQ(q => q.root.subtree().withTag('branched-snfg-3d'));
    for (const node of snfgNodes) {
        if (!node.state.isHidden) {
            setSubtreeVisibility(plugin.state.data, node.transform.ref, true); // true means hidden
            hiddenNodes.set(node.transform.ref, true);
        }
    }
}

/** Make visible all carbohydrate SNFG visuals that have been hidden by `hideSnfgNodes` */
function unhideSnfgNodes(plugin: PluginContext) {
    const hiddenNodes = getExtensionState(plugin).hiddenNodes;
    for (const nodeRef of Array.from(hiddenNodes.keys())) {
        try {
            setSubtreeVisibility(plugin.state.data, nodeRef, false); // false means visible
        } catch {
            // this is OK, the node has been removed
        }
    }
    hiddenNodes.clear();
}


interface ExtensionState {
    structConnCache: NodeMap<{ [id: string]: StructConnRecord }>,
    hiddenNodes: NodeMap<boolean>,
}
function createExtensionState(plugin: PluginContext): ExtensionState {
    return {
        structConnCache: new NodeMap<{ [id: string]: StructConnRecord }>(plugin),
        hiddenNodes: new NodeMap<boolean>(plugin),
    };
}
function getExtensionState(plugin: PluginContext): ExtensionState {
    return (plugin.customState as any).StructConn ??= createExtensionState(plugin);
}


/** Map in which keys are refs of plugin state tree nodes,
 * and they get automatically deleted from the map
 * when the node is deleted from the tree. */
class NodeMap<T> extends PluginBehavior.WithSubscribers<undefined> {
    private readonly map = new Map<string, T>();
    private sub?: PluginCommand.Subscription = undefined;

    constructor(plugin: PluginContext) {
        super(plugin, undefined);
    }
    register(): void { }
    keys() {
        return this.map.keys();
    }
    has(key: string) {
        return this.map.has(key);
    }
    get(key: string) {
        return this.map.get(key);
    }
    set(key: string, value: T) {
        this.sub ??= this.subscribeObservable(this.plugin.state.data.events.cell.removed, o => this.delete(o.ref));
        this.map.set(key, value);
    }
    delete(key: string) {
        const deleted = this.map.delete(key);
        if (this.map.size === 0 && this.sub) {
            this.sub.unsubscribe();
            this.sub = undefined;
        }
        return deleted;
    }
    clear() {
        this.map.clear();
        if (this.sub) {
            this.sub.unsubscribe();
            this.sub = undefined;
        }
    }
}
