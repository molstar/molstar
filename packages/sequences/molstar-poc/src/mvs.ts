/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Structures shown through MolViewSpec (seqspec:libraries#the-mol-adapter): what Seq* asks the 3D view for (this
// structure, these chains, these colors, this selection) is an MVS state, never ad-hoc Mol* calls. Colors travel as MVS
// annotation tables (the `residue_range` schema, columns as arrays: label_asym_id, beg/end_label_seq_id, color), each
// applied by a `color_from_uri` node over a `data:` URI, in order, so a later table wins where two color one residue.

import type { MVSData_State } from '@molstar/mvs-builder/mvs-data';
import { createMVSBuilder } from '@molstar/mvs-builder/tree/mvs/mvs-builder';
import type { ColorT } from '@molstar/mvs-builder/tree/mvs/param-types';

/** An MVS state (Mol*'s `MVSData_State`), what a host loads. */
export type MvsState = MVSData_State;

/** A chain by MVS's selector fields; without `instance_id`, every copy. */
export interface MvsChain {
  readonly label_asym_id: string;
  readonly instance_id?: string;
}
/** An MVS annotation in the `residue_range` schema, object of arrays (MVS's compact form). */
export interface MvsColorTable {
  readonly label_asym_id: string[];
  readonly beg_label_seq_id: number[];
  readonly end_label_seq_id: number[];
  readonly color: string[];
}
export const MvsColorTable = {
  create: (): MvsColorTable => ({ label_asym_id: [], beg_label_seq_id: [], end_label_seq_id: [], color: [] }),
  size: (t: MvsColorTable) => t.color.length,
  /** As a `data:` URI, what `color_from_uri` reads (a demonstration; a host may serve the JSON instead). */
  dataUri: (t: MvsColorTable) => `data:application/json;base64,${base64(JSON.stringify(t))}`,
};
function base64(text: string): string {
  const bytes = new TextEncoder().encode(text);
  let s = '';
  for (let i = 0; i < bytes.length; i++) s += String.fromCharCode(bytes[i]);
  return btoa(s);
}

/**
 * What a structure's own defaults say (`custom.molstar` on an annotation, seqspec:format#extensions): MVS template
 * values.
 */
export interface MolstarDefaults {
  readonly representation?: 'cartoon' | 'surface' | 'ball_and_stick' | 'spacefill' | 'backbone';
  /** The polymer's color where no table reaches. */
  readonly base?: string;
  /** Ligands and ions as ball and stick (default true). */
  readonly ligands?: boolean;
  /** The other chains' opacity when some are in front (default 0.25). */
  readonly faint?: number;
}

export interface StructureScene {
  readonly url: string;
  readonly format: 'bcif' | 'mmcif' | 'pdb';
  /** The model, 0-based (`model_index`). */
  readonly model?: number;
  /** An assembly id; null: the deposited model. */
  readonly assembly?: string | null;
  /** The chains shown in front, the others faint; null: all in front. */
  readonly chains?: readonly MvsChain[] | null;
  /** Color tables in order: a later one wins. */
  readonly colors?: readonly MvsColorTable[];
  /** Residues shown as ball and stick too, each table's in its colors (a site track's: PTMs, binding sites). */
  readonly sticks?: readonly MvsColorTable[];
  /** Residues shown as ball and stick too, with a label (the Seq* selection). */
  readonly selection?: { readonly residues: MvsColorTable; readonly label: string } | null;
  readonly defaults?: MolstarDefaults;
  /** The camera on `chains` (default true when chains are given). */
  readonly focus?: boolean;
  /** The canvas's background (MVS's `canvas` node; a host's light or dark surface); unset: MVS's default. */
  readonly background?: string;
}

const rows = (t: MvsColorTable) =>
  t.color.map((_, i) => ({
    label_asym_id: t.label_asym_id[i],
    beg_label_seq_id: t.beg_label_seq_id[i],
    end_label_seq_id: t.end_label_seq_id[i],
  }));

/** One MVS state: the structure in its representation (its chains in front), colored by the tables, the selection. */
export function structureMvs(scene: StructureScene): MVSData_State {
  const d = scene.defaults ?? {};
  const b = createMVSBuilder();
  if (scene.background) b.canvas({ background_color: scene.background as ColorT });
  const parsed = b.download({ url: scene.url }).parse({ format: scene.format });
  const s =
    scene.assembly != null
      ? parsed.assemblyStructure({ assembly_id: scene.assembly, model_index: scene.model ?? 0 })
      : parsed.modelStructure({ model_index: scene.model ?? 0 });
  const chains = scene.chains?.length ? scene.chains.map((c) => ({ ...c })) : null;
  const front = s.component({ selector: chains ?? 'polymer' });
  const repr = front.representation({ type: d.representation ?? 'cartoon' });
  repr.color({ color: (d.base ?? '#c8c8c8') as ColorT });
  for (const t of scene.colors ?? []) {
    if (MvsColorTable.size(t) === 0) continue;
    repr.colorFromUri({ uri: MvsColorTable.dataUri(t), format: 'json', schema: 'residue_range' });
  }
  if (chains) {
    s.component({ selector: 'polymer' })
      .representation({ type: d.representation ?? 'cartoon' })
      .color({ color: '#e0e0e0' })
      .opacity({ opacity: d.faint ?? 0.25 });
    if (scene.focus !== false) front.focus({});
  }
  if (d.ligands !== false) {
    s.component({ selector: 'ligand' }).representation({ type: 'ball_and_stick' }).color({ color: '#7b6f62' });
    s.component({ selector: 'ion' }).representation({ type: 'ball_and_stick' });
  }
  for (const t of scene.sticks ?? []) {
    if (MvsColorTable.size(t) === 0) continue;
    const uri = MvsColorTable.dataUri(t);
    s.componentFromUri({ uri, format: 'json', schema: 'residue_range' })
      .representation({ type: 'ball_and_stick' })
      .colorFromUri({ uri, format: 'json', schema: 'residue_range' });
  }
  const sel = scene.selection;
  if (sel && MvsColorTable.size(sel.residues) > 0) {
    const c = s.component({ selector: rows(sel.residues) });
    c.representation({ type: 'ball_and_stick' }).color({ color: (sel.residues.color[0] ?? '#e8a33d') as ColorT });
    c.label({ text: sel.label });
  }
  return b.getState();
}
