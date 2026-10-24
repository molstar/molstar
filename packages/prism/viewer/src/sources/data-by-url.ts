/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CifFrame } from '@molstar/io/reader/cif';
import type { Trajectory } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
import { PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { Download, ReadFile } from '@molstar/plugin/state/transforms/data/fetch';
import type { StateObjectCell } from '@molstar/core/state';
import { Asset } from '@molstar/core/util/assets';
import { StringLike } from '@molstar/core/util/string-like';
import { merge } from 'rxjs';
import type { DataByUrl, DataByUrlEntry } from '../sequence/model.js';
import type { StructureFormat, StructureSource } from '../sequence/molstar/source.js';
import type { SourceScope } from './state.js';

// The sequence's parse cache (sequence.md §2): `dataByUrl`, the plugin's data cells by address (a Download's URL, a
// read file's name) with the parse under each, so Seq*'s loaders reuse what the plugin holds instead of refusing or
// reading it again. Kept lazily: a data-tree change only marks it stale (O(1) per event), the next lookup rebuilds it
// (O(cells), edge rate: a document's load), so a scene of thousands of cells pays nothing while no document loads.

/** The structure formats the sequence's adapter reads, by the trajectory transformer that parsed them. */
const Trajectories: Readonly<Record<string, StructureFormat | 'cif'>> = {
  'trajectory-from-mmcif': 'cif',
  'trajectory-from-pdb': 'pdb',
  'trajectory-from-sdf': 'sdf',
  'trajectory-from-mol2': 'mol2',
};

/** How deep under a data cell its parse is looked for (Download → ParseCif → TrajectoryFromMmCif). */
const ParseDepth = 3;

interface Indexed extends DataByUrlEntry {
  readonly cell: StateObjectCell;
  readonly trajectoryCell?: StateObjectCell;
}

export function createDataByUrl(scope: SourceScope, plugin: PluginContext): DataByUrl {
  let index: Map<string, Indexed> | undefined;
  /** Sources made per trajectory: one `StructureSource` for every request of it (its structures made once). */
  const sources = new WeakMap<Trajectory, StructureSource>();
  const events = plugin.state.data.events;
  scope.add(
    merge(events.object.created, events.object.updated, events.object.removed, events.cell.removed).subscribe(() => {
      index = undefined;
    }),
  );

  const parseUnder = (ref: string): StateObjectCell | undefined => {
    const { tree, cells } = plugin.state.data;
    let level = [ref];
    for (let d = 0; d < ParseDepth && level.length; d++) {
      const next: string[] = [];
      for (const r of level) {
        tree.children.get(r)?.forEach((c) => {
          next.push(c);
        });
      }
      for (const r of next) {
        const cell = cells.get(r);
        const name = cell?.transform.transformer.definition.name;
        if (cell?.obj && name && Trajectories[name] && SO.Molecule.Trajectory.is(cell.obj)) return cell;
      }
      level = next;
    }
    return undefined;
  };

  const build = () => {
    const map = new Map<string, Indexed>();
    for (const cell of plugin.state.data.cells.values()) {
      const t = cell.transform.transformer;
      const params = cell.transform.params as { url?: string | Asset.Url; file?: Asset.File | null } | undefined;
      const address =
        t === Download && params?.url
          ? Asset.getUrl(params.url)
          : t === ReadFile && params?.file
            ? params.file.name
            : undefined;
      const data = cell.obj?.data;
      const held = typeof data === 'string' || data instanceof Uint8Array || data instanceof SO.Data.ResponseBody;
      if (!address || !held || map.has(address)) continue;
      const trajectoryCell = parseUnder(cell.transform.ref);
      map.set(address, {
        address,
        data: cell.transform.ref,
        ...(trajectoryCell ? { parse: trajectoryCell.transform.ref, trajectoryCell } : {}),
        cell,
      });
    }
    return map;
  };
  const lookup = (address: string) => (index ??= build()).get(address);

  return {
    entry(address) {
      const e = lookup(address);
      return e && { address: e.address, data: e.data, ...(e.parse ? { parse: e.parse } : {}) };
    },
    bytes(address) {
      const data = lookup(address)?.cell.obj?.data;
      if (data instanceof Uint8Array) return data;
      // a response's held bytes, as a binary cell's (no copy); one holding text is encoded as a string cell is, and
      // keeps its form
      if (data instanceof SO.Data.ResponseBody) {
        return data.heldBytes ?? new TextEncoder().encode(StringLike.toString(data.textNow()));
      }
      return typeof data === 'string' ? new TextEncoder().encode(data) : undefined;
    },
    structure(address) {
      const e = lookup(address);
      const tc = e?.trajectoryCell;
      const trajectory = tc?.obj?.data as Trajectory | undefined;
      if (!e || !tc || !trajectory) return undefined;
      const hit = sources.get(trajectory);
      if (hit) return hit;
      const kind = Trajectories[tc.transform.transformer.definition.name];
      const binary = e.cell.obj?.data instanceof Uint8Array;
      const format: StructureFormat = kind === 'cif' ? (binary ? 'bcif' : 'mmcif') : kind;
      const frame = cifFrameOf(plugin, tc);
      const source: StructureSource = {
        kind: 'molstar:structure-source',
        format,
        trajectory,
        frame,
        structures: new Map(),
      };
      sources.set(trajectory, source);
      return source;
    },
  };
}

/**
 * The CIF data block a trajectory was parsed from (the parse cell above it, the block its params name), for the
 * categories the sequence reads beside the structure: `ma_qa_metric_local` (a predicted model), SIFTS numberings,
 * `pdbx_poly_seq_scheme` (sequence.md §2). Null: not from a CIF, or every block in one trajectory.
 */
export function cifFrameOf(plugin: PluginContext, trajectoryCell: StateObjectCell): CifFrame | null {
  const parent = plugin.state.data.cells.get(trajectoryCell.transform.parent);
  const cif = parent?.obj && SO.Format.Cif.is(parent.obj) ? parent.obj.data : undefined;
  if (!cif?.blocks.length) return null;
  const p = trajectoryCell.transform.params as
    | { loadAllBlocks?: boolean; blockHeader?: string; blockIndex?: number }
    | undefined;
  if (p?.loadAllBlocks) return null;
  const header = p?.blockHeader?.toLowerCase();
  const block = header
    ? cif.blocks.find((b) => b.header.toLowerCase() === header)
    : cif.blocks[Math.max(0, p?.blockIndex ?? 0)];
  return (block as CifFrame | undefined) ?? null;
}
