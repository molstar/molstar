/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateSelection } from '@molstar/core/state';
import { Task } from '@molstar/core/task';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Structure } from '@molstar/model/model/structure';
import type { CustomProperty } from '@molstar/model/props/common/custom-property';
import type { PluginContext } from '@molstar/plugin/context';
import type { StructureRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import { type ApiContext, ApiArea, type ApiFailure, ApiResult } from '@molstar/prism-viewer/api/area';
import { StructureParam } from '@molstar/prism-viewer/api/params';
import { isFailure, labelOf, oneStructure } from '@molstar/prism-viewer/api/scope';
import { watchedState } from '@molstar/prism-viewer/sources/watched';
import { type ObservableState, observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import {
  AssemblySymmetryData,
  AssemblySymmetryDataProvider,
  type AssemblySymmetryDataValue,
  AssemblySymmetryProvider,
} from './prop.js';
import { AssemblySymmetry3D, getAssemblySymmetryConfig } from './transforms.js';

// The `assembly-symmetry` area (api.md §2, parity.md §6): RCSB PDB's annotation of an assembly's symmetry
// (BioJava's, from data.rcsb.org; PDBe's when the viewer is configured so), its axes in 3D and the coloring by its
// clusters. `get` lists the symmetries, `show` adds the axes for one (the first that is not C1 by default), `hide`
// removes them, `color` colors the structure by the symmetry's clusters. Each reaches only the configured service, so
// `addresses` is empty (agent.md §4). The chosen symmetry is kept in the structure's custom-properties node, as v5's
// controls did, so Undo and sessions carry it.

const Tag = AssemblySymmetryData.Tag;
const PropName = AssemblySymmetryProvider.descriptor.name;

/** One symmetry of an assembly, as `get` lists it. */
export interface SymmetryInfo {
  /** Its number, from 1: what `show` and `color` take. */
  readonly index: number;
  /** "Global Symmetry", "Pseudo Symmetry", "Local Symmetry". */
  readonly kind: string;
  /** The point group: C1, C2, D3, H, … */
  readonly symbol: string;
  /** "Cyclic", "Dihedral", "Helical", … */
  readonly type: string;
  readonly oligomericState: string;
  readonly stoichiometry: readonly string[];
  readonly clusters: number;
  readonly axes: number;
}

export interface AssemblySymmetryInfo {
  readonly structure: string;
  readonly symmetries: readonly SymmetryInfo[];
  /** The symmetry whose axes are shown, by number. */
  readonly shown?: number;
}

const network = { modules: ['network'] as const, approval: 'network' as const, addresses: () => [] };

const SymmetryParam = () =>
  PD.Numeric(
    0,
    { min: 0, step: 1 },
    { description: 'A symmetry by its number from assembly-symmetry.get (1, 2, …); 0: the first that is not C1' },
  );

/** What the assembly symmetry's row shows of a structure (`prism-viewer-ui.tsx`): read, never fetched for it. */
export interface SymmetryRowState {
  /** A PDB entry's assembly: RCSB PDB annotates its symmetry. */
  readonly applicable: boolean;
  /** The symmetries, once asked for (`get`, `show` or `color` fetched them). */
  readonly symmetries?: readonly SymmetryInfo[];
  /** The symmetry chosen (by number), its axes shown or not. */
  readonly selected?: number;
  /** Its axes are in the scene. */
  readonly shown: boolean;
}

export class AssemblySymmetryApi extends ApiArea {
  readonly id = 'assembly-symmetry';

  /**
   * The row's state for a structure (by state ref), following the state tree while watched (a view bound to it): a
   * read for views, as the animation area's `state` (the Prism viewer's `_spec/README.md`, Folders). It never fetches:
   * the symmetries show once an operation has asked RCSB PDB for them.
   */
  /**
   * Whether a structure (by state ref) is a PDB entry's assembly, which RCSB PDB annotates: the row's `when`, so the
   * Structure panel gives a row only to the structures it is for (its rows are of one height: a row that renders
   * nothing would leave a gap).
   */
  applies(structure: string): boolean {
    return isApplicableRef(this.ctx.plugin, structure);
  }

  rowState(structure: string): ReadonlyState<SymmetryRowState> {
    const plugin = this.ctx.plugin;
    return watchedState(
      () => readRow(plugin, structure),
      // A fetch changes no state (a C1-only answer adds no axes): the row hears of it apart.
      () => [plugin.state.data.events.changed, fetchedSignal(plugin).asReadonly],
      { compare: (a, b) => JSON.stringify(a) === JSON.stringify(b) },
    );
  }

  readonly get = this.operation({
    id: 'assembly-symmetry.get',
    title: 'Describe the assembly symmetry',
    params: { structure: StructureParam() },
    read: true,
    ...network,
    async run({ structure }, ctx) {
      const s = applicable(ctx, structure);
      if (isFailure(s)) return s;
      const list = await fetchList(ctx, s);
      if (isFailure(list)) return list;
      const info = describe(ctx, s, list);
      return ApiResult.ok(summary(info), {}, { value: info });
    },
  });

  readonly show = this.operation({
    id: 'assembly-symmetry.show',
    title: 'Show the symmetry axes',
    params: { symmetry: SymmetryParam(), structure: StructureParam() },
    ...network,
    async run({ symmetry, structure }, ctx) {
      const s = applicable(ctx, structure);
      if (isFailure(s)) return s;
      const picked = await pick(ctx, s, symmetry);
      if (isFailure(picked)) return picked;
      const name = labelOf(ctx, s);
      const said = `${picked.kind}: ${picked.symbol} (${picked.type.toLowerCase()})`;
      const axes = axesRef(ctx.plugin, s);
      if (axes && selected(s) === picked.index - 1) {
        return ApiResult.already(`${name} already shows the axes of its ${said}.`, { value: picked });
      }
      const failed = await apply(ctx, s, picked.index - 1, true);
      if (failed) return failed;
      return ApiResult.ok(
        `Showed the axes of ${name}'s ${said}.`,
        { nodes: [s.cell.transform.ref] },
        { value: picked },
      );
    },
  });

  readonly hide = this.operation({
    id: 'assembly-symmetry.hide',
    title: 'Hide the symmetry axes',
    params: { structure: StructureParam() },
    async run({ structure }, ctx) {
      const s = oneStructure(ctx, structure || undefined);
      if (isFailure(s)) return s;
      const ref = axesRef(ctx.plugin, s);
      if (!ref) return ApiResult.already(`${labelOf(ctx, s)} shows no symmetry axes.`);
      await ctx.plugin.state.data.build().delete(ref).commit();
      return ApiResult.ok(`Hid the symmetry axes of ${labelOf(ctx, s)}.`, { nodes: [s.cell.transform.ref] });
    },
  });

  readonly color = this.operation({
    id: 'assembly-symmetry.color',
    title: 'Color by symmetry cluster',
    params: { symmetry: SymmetryParam(), structure: StructureParam() },
    ...network,
    async run({ symmetry, structure }, ctx) {
      const s = applicable(ctx, structure);
      if (isFailure(s)) return s;
      const picked = await pick(ctx, s, symmetry);
      if (isFailure(picked)) return picked;
      if (!picked.clusters) {
        return ApiResult.fail('blocked', `The ${picked.symbol} symmetry of ${labelOf(ctx, s)} has no clusters.`);
      }
      if (selected(s) !== picked.index - 1) {
        const failed = await apply(ctx, s, picked.index - 1, false);
        if (failed) return failed;
      }
      // The theme reads the chosen symmetry's clusters (`color.ts`); the coloring itself is the core's.
      const r = await ctx.api.representation.color({
        target: 'all',
        theme: Tag.Cluster,
        structure: s.cell.transform.ref,
      });
      if (!r.ok) return r;
      return ApiResult.ok(`Colored ${labelOf(ctx, s)} by the clusters of its ${picked.symbol} symmetry.`, r.changed, {
        value: picked,
      });
    },
  });
}

// --- helpers ----------------------------------------------------------------------------------------------------

function isApplicableRef(plugin: PluginContext, ref: string): boolean {
  const s = plugin.managers.structure.hierarchy.current.refs.get(ref);
  const data = s?.kind === 'structure' ? s.cell.obj?.data : undefined;
  return !!data && AssemblySymmetryData.isApplicable(data);
}

function readRow(plugin: PluginContext, ref: string): SymmetryRowState {
  const s = plugin.managers.structure.hierarchy.current.refs.get(ref);
  const data = s?.kind === 'structure' ? s.cell.obj?.data : undefined;
  if (s?.kind !== 'structure' || !data || !AssemblySymmetryData.isApplicable(data)) {
    return { applicable: false, shown: false };
  }
  const list = AssemblySymmetryDataProvider.get(data).value;
  const index = selected(s);
  return {
    applicable: true,
    shown: !!axesRef(plugin, s),
    ...(list ? { symmetries: infos(list) } : {}),
    ...(index >= 0 ? { selected: index + 1 } : {}),
  };
}

function infos(list: AssemblySymmetryDataValue): SymmetryInfo[] {
  return list.map((x, i) => ({
    index: i + 1,
    kind: x.kind,
    symbol: x.symbol,
    type: x.type,
    oligomericState: x.oligomeric_state,
    stoichiometry: x.stoichiometry,
    clusters: x.clusters?.length ?? 0,
    axes: x.rotation_axes?.length ?? 0,
  }));
}

function applicable(ctx: ApiContext, name: string): StructureRef | ApiFailure {
  const s = oneStructure(ctx, name || undefined);
  if (isFailure(s)) return s;
  const data = s.cell.obj?.data;
  if (!data || !AssemblySymmetryData.isApplicable(data)) {
    return ApiResult.fail(
      'blocked',
      `RCSB PDB annotates the symmetry of a PDB entry's assemblies; ${labelOf(ctx, s)} is not one.`,
      'Load the entry from the PDB (load.entry) and show one of its assemblies (structure.assembly).',
    );
  }
  return s;
}

function propertyContext(plugin: PluginContext, runtime: CustomProperty.Context['runtime']): CustomProperty.Context {
  return { runtime, assetManager: plugin.managers.asset, errorContext: plugin.errorContext };
}

const fetches = new WeakMap<PluginContext, ObservableState<number>>();

/** Bumped after each fetch, answered or failed: the rows' states read again. */
function fetchedSignal(plugin: PluginContext): ObservableState<number> {
  let s = fetches.get(plugin);
  if (!s) fetches.set(plugin, (s = observableState(0)));
  return s;
}

/** The symmetries the service has for the structure's assembly (fetched once per structure). */
async function fetchList(ctx: ApiContext, s: StructureRef): Promise<AssemblySymmetryDataValue | ApiFailure> {
  const data = s.cell.obj!.data;
  const { DefaultServerType: serverType, DefaultServerUrl: serverUrl } = getAssemblySymmetryConfig(ctx.plugin);
  try {
    await ctx.plugin.runTask(
      Task.create('Assembly symmetry', (runtime) =>
        AssemblySymmetryDataProvider.attach(propertyContext(ctx.plugin, runtime), data, { serverType, serverUrl }),
      ),
    );
  } catch (e) {
    fetchedSignal(ctx.plugin).set(fetchedSignal(ctx.plugin).value + 1);
    return ApiResult.fail(
      'failed',
      `Could not get the symmetry of ${labelOf(ctx, s)} from ${serverType === 'pdbe' ? 'PDBe' : 'RCSB PDB'}: ${message(e)}`,
      'Try again later.',
    );
  }
  fetchedSignal(ctx.plugin).set(fetchedSignal(ctx.plugin).value + 1);
  return AssemblySymmetryDataProvider.get(data).value ?? [];
}

function describe(ctx: ApiContext, s: StructureRef, list: AssemblySymmetryDataValue): AssemblySymmetryInfo {
  const shown = axesRef(ctx.plugin, s) ? selected(s) : -1;
  return {
    structure: labelOf(ctx, s),
    symmetries: infos(list),
    ...(shown >= 0 ? { shown: shown + 1 } : {}),
  };
}

function summary(info: AssemblySymmetryInfo): string {
  if (!info.symmetries.length) return `RCSB PDB has no symmetry for ${info.structure}.`;
  const lines = info.symmetries.map(
    (x) =>
      `${x.index}. ${x.kind}: ${x.symbol} (${x.type.toLowerCase()}), ${x.oligomericState}, stoichiometry ${x.stoichiometry.join(', ')}${info.shown === x.index ? '; its axes are shown' : ''}`,
  );
  return [`${info.structure}:`, ...lines].join('\n');
}

/** The symmetry asked for (by number; 0: the first that is not C1), with what the service has. */
async function pick(ctx: ApiContext, s: StructureRef, n: number): Promise<SymmetryInfo | ApiFailure> {
  const list = await fetchList(ctx, s);
  if (isFailure(list)) return list;
  const info = describe(ctx, s, list);
  const name = labelOf(ctx, s);
  if (!n) {
    const first = info.symmetries.find((x) => x.symbol !== 'C1');
    if (first) return first;
    return ApiResult.fail(
      'blocked',
      info.symmetries.length
        ? `${name} has no symmetry but C1: there are no axes.`
        : `RCSB PDB has no symmetry for ${name}.`,
    );
  }
  const chosen = info.symmetries[n - 1];
  if (!chosen) {
    return ApiResult.fail(
      'invalid',
      `${name} has ${info.symmetries.length} symmetr${info.symmetries.length === 1 ? 'y' : 'ies'}.`,
      `Give 1 to ${info.symmetries.length}, or 0 for the first that is not C1.`,
    );
  }
  if (chosen.symbol === 'C1') {
    return ApiResult.fail('invalid', `Symmetry ${n} of ${name} is C1: it has no axes or clusters.`);
  }
  return chosen;
}

/** The symmetry the structure has chosen (the property's index; -1: none). */
function selected(s: StructureRef): number {
  const data = s.cell.obj?.data as Structure | undefined;
  if (!data || !AssemblySymmetryProvider.get(data.root).value) return -1;
  return AssemblySymmetryProvider.props(data.root).symmetryIndex;
}

function axesRef(plugin: PluginContext, s: StructureRef): string | undefined {
  return StateSelection.findTagInSubtree(plugin.state.data.tree, s.cell.transform.ref, Tag.Representation);
}

/**
 * Chooses the symmetry, in the structure's custom-properties node where it has one (so the state holds it), and adds
 * or updates the axes when asked; one state update.
 */
async function apply(ctx: ApiContext, s: StructureRef, index: number, axes: boolean): Promise<ApiFailure | undefined> {
  const plugin = ctx.plugin;
  const { DefaultServerType: serverType, DefaultServerUrl: serverUrl } = getAssemblySymmetryConfig(plugin);
  const values = { serverType, serverUrl, symmetryIndex: index };
  const data = s.cell.obj!.data;
  const props = s.properties?.cell;
  const build = plugin.state.data.build();
  if (props) {
    const p = props.transform.params as { autoAttach: string[]; properties: Record<string, unknown> };
    build.to(props.transform.ref).update({
      autoAttach: p.autoAttach.includes(PropName) ? p.autoAttach : [...p.autoAttach, PropName],
      properties: { ...p.properties, [PropName]: values },
    });
  } else {
    await plugin.runTask(
      Task.create('Assembly symmetry', (runtime) =>
        AssemblySymmetryProvider.attach(propertyContext(plugin, runtime), data, values),
      ),
    );
    // Without the node, the axes read the property when they are made: made again for another symmetry.
    const old = axesRef(plugin, s);
    if (old && axes) build.delete(old);
  }
  if (axes) build.to(s.cell.transform.ref).applyOrUpdateTagged(Tag.Representation, AssemblySymmetry3D);
  try {
    await build.commit({ revertOnError: true });
  } catch (e) {
    return ApiResult.fail('failed', `Could not show the symmetry of ${labelOf(ctx, s)}: ${message(e)}`);
  }
  if (selected(s) !== index) {
    return ApiResult.fail('failed', `Could not get symmetry ${index + 1} of ${labelOf(ctx, s)}.`, 'Try again later.');
  }
  return undefined;
}

const message = (e: unknown) => (e instanceof Error ? e.message : String(e));
