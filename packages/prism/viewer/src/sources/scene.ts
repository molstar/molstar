/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Volume } from '@molstar/model/model/volume';
import type { ParticleList } from '@molstar/model/model/particles/particle-list';
import {
  type Structure,
  StructureElement,
  StructureProperties,
  StructureSelection,
  type Unit,
} from '@molstar/model/model/structure';
import { Script } from '@molstar/model/script/script';
import { MolScriptBuilder as MS } from '@molstar/query-language/language/builder';
import type { PluginContext } from '@molstar/plugin/context';
import type {
  StructureComponentRef,
  StructureRef,
  StructureRepresentationRef,
} from '@molstar/plugin/state/manager/structure/hierarchy-state';
import type {
  ParticleListRef,
  ParticleRepresentationRef,
} from '@molstar/plugin/state/manager/particles/hierarchy-state';
import type { VolumeRef, VolumeRepresentationRef } from '@molstar/plugin/state/manager/volume/hierarchy-state';
import { type IconName, type ObservableState, observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import type { ViewerApi } from '../api/area.js';
import type { SceneInspection, SceneMenuItem, SceneModel, SceneNode, SceneNodeKind } from '../model/scene.js';
import { namesOf } from '../api/scope.js';
import { gridText, shortNumber, voxelText, volumeStats } from '../api/volume.js';
import { structureNames } from './naming.js';
import type { SourceScope } from './state.js';

// The scene as the Structure tab shows it (state.md §1, row 2; view-models.md §3): one record per node, keyed by state
// ref, so a node keeps its identity across hierarchy rebuilds, and each field its own state, so a row updates only where
// it changed. A structure's contents (entities, then chains) are computed when first expanded. Every action goes through
// the API, bound to the person (README §3.3): the tree changes state only by operations.

const ReprIcons: Readonly<Record<string, IconName>> = {
  cartoon: 'cartoon',
  'ball-and-stick': 'ball-and-stick',
  spacefill: 'spacefill',
  'molecular-surface': 'surface',
  'gaussian-surface': 'surface',
  'gaussian-volume': 'surface',
  putty: 'cartoon',
  backbone: 'cartoon',
  line: 'ball-and-stick',
  label: 'tag',
  carbohydrate: 'shapes',
  ellipsoid: 'spacefill',
};

/** The operation that applies an inspector form, by the kind of representation it edits. */
const UpdateOperations: Partial<Record<SceneNodeKind, string>> = {
  representation: 'representation.update',
  'volume-representation': 'volume.update',
  'particles-representation': 'particles.update',
};

class Node implements SceneNode {
  readonly #label: ObservableState<string>;
  readonly #description = observableState<string | undefined>(undefined);
  readonly #hidden = observableState(false);
  readonly #expanded: ObservableState<boolean>;
  readonly #children = observableState<readonly SceneNode[]>([]);
  icon: IconName;
  menu?: readonly SceneMenuItem[];
  removable?: boolean;
  lazy?: boolean;
  load?: () => void;

  constructor(
    readonly id: string,
    readonly kind: SceneNodeKind,
    icon: IconName,
    label: string,
    expanded: boolean,
    readonly actions: {
      visibility?(hidden: boolean): void;
      highlight?(on: boolean): void;
      focus?(): void;
      remove?(): void;
    },
  ) {
    this.icon = icon;
    this.#label = observableState(label);
    this.#expanded = observableState(expanded);
  }

  get label(): ReadonlyState<string> {
    return this.#label.asReadonly;
  }
  get description(): ReadonlyState<string | undefined> {
    return this.#description.asReadonly;
  }
  get isHidden(): ReadonlyState<boolean> {
    return this.#hidden.asReadonly;
  }
  get isExpanded(): ReadonlyState<boolean> {
    return this.#expanded.asReadonly;
  }
  get children(): ReadonlyState<readonly SceneNode[]> {
    return this.#children.asReadonly;
  }

  set(label: string, description: string | undefined, hidden: boolean): void {
    this.#label.set(label);
    this.#description.set(description);
    this.#hidden.set(hidden);
  }

  setChildren(children: readonly SceneNode[]): void {
    this.#children.set(children);
  }

  toggleVisibility(): void {
    this.actions.visibility?.(!this.#hidden.value);
  }
  toggleExpanded(): void {
    this.#expanded.set(!this.#expanded.value);
    if (this.#expanded.value && this.lazy) this.load?.();
  }
  highlight(on: boolean): void {
    this.actions.highlight?.(on);
  }
  focus(): void {
    this.actions.focus?.();
  }
  remove(): void {
    this.actions.remove?.();
  }
}

export interface SceneSource extends SceneModel {
  /** Rebuilds from the hierarchy (called on its changes; tests may call it). */
  sync(): void;
}

/**
 * The scene over the plugin. `api` is bound to the person (its calls are the person's steps); `scope` holds the
 * subscriptions.
 */
export function createSceneSource(scope: SourceScope, plugin: PluginContext, api: ViewerApi): SceneSource {
  const nodes = new Map<string, Node>();
  const roots = observableState<readonly SceneNode[]>([]);
  const picked = observableState<SceneNode | undefined>(undefined);
  const inspections = new Map<string, SceneInspection>();

  const structureNode = (s: StructureRef, name: string): Node => {
    const ref = s.cell.transform.ref;
    let n = nodes.get(ref);
    if (!n) {
      n = new Node(ref, 'structure', 'structure', name, true, {
        visibility: (hide) => void api.run('structure.visibility', { structure: ref, visible: !hide }),
        highlight: (on) =>
          on ? void api.run('highlight.set', { target: 'all', structure: ref }) : void api.run('highlight.clear'),
        focus: () => void api.run('focus.set', { target: 'all', structure: ref }),
        remove: () => void api.run('structure.remove', { structure: ref }),
      });
      n.removable = true;
      nodes.set(ref, n);
    }
    const data = s.cell.obj?.data;
    n.set(name, data ? describeStructure(data) : undefined, !!s.cell.state.isHidden);
    return n;
  };

  const componentNode = (c: StructureComponentRef, structureRef: string): Node => {
    const ref = c.cell.transform.ref;
    let n = nodes.get(ref);
    if (!n) {
      const loci = () => {
        const data = c.cell.obj?.data;
        return data ? Structure_toLoci(data) : undefined;
      };
      n = new Node(ref, 'component', 'layers', c.cell.obj?.label ?? 'Component', false, {
        visibility: (hide) =>
          void api.run('component.visibility', { component: ref, structure: structureRef, visible: !hide }),
        highlight: (on) => {
          const l = loci();
          if (on && l) void api.run('highlight.set', { target: { loci: l } });
          else void api.run('highlight.clear');
        },
        focus: () => {
          const l = loci();
          if (l) void api.run('focus.set', { target: { loci: l } });
        },
        remove: () => void api.run('component.remove', { component: ref, structure: structureRef }),
      });
      n.removable = true;
      n.menu = componentMenu(plugin, api, ref, structureRef, loci);
      nodes.set(ref, n);
    }
    const reprs = c.representations.map((r) => reprLabel(plugin, r)).join(', ');
    const atoms = c.cell.obj?.data.elementCount;
    n.set(
      c.cell.obj?.label ?? 'Component',
      [reprs, atoms !== undefined ? `${atoms.toLocaleString('en')} atoms` : undefined].filter(Boolean).join(' · ') ||
        undefined,
      !!c.cell.state.isHidden,
    );
    return n;
  };

  const reprNode = (r: StructureRepresentationRef, componentRef: string, structureRef: string): Node => {
    const ref = r.cell.transform.ref;
    let n = nodes.get(ref);
    const type = (r.cell.transform.params as { type?: { name: string } } | undefined)?.type?.name ?? '';
    if (!n) {
      n = new Node(ref, 'representation', ReprIcons[type] ?? 'cartoon', reprLabel(plugin, r), false, {
        visibility: (hide) => void api.run('representation.visibility', { representation: ref, visible: !hide }),
        remove: () => void api.run('representation.remove', { component: componentRef, type, structure: structureRef }),
      });
      n.removable = true;
      nodes.set(ref, n);
    }
    n.icon = ReprIcons[type] ?? 'cartoon';
    n.set(reprLabel(plugin, r), colorLabel(plugin, r), !!r.cell.state.isHidden);
    // An open inspector follows the cell (Undo, another surface), without reporting it as an edit.
    inspections.get(ref)?.form?.load(r.cell.transform.params as never);
    return n;
  };

  /** A structure's contents: entities, then their chains; computed when first expanded. */
  const contentsNode = (s: StructureRef, structureRef: string): Node => {
    const id = `${structureRef}#contents`;
    let n = nodes.get(id);
    if (!n) {
      const node = new Node(id, 'contents', 'list-tree', 'Contents', false, {});
      node.lazy = true;
      node.removable = false;
      node.load = () => {
        node.lazy = false;
        const data = s.cell.obj?.data;
        if (!data) return;
        node.setChildren(entities(data, structureRef, nodes, api));
      };
      n = node;
      nodes.set(id, n);
    }
    const data = s.cell.obj?.data;
    n.set('Contents', data ? `${data.model.entities.data._rowCount} entities` : undefined, false);
    return n;
  };

  // Volumes and particles (plan 2.5): top-level nodes like structures, their representations as children. Eyes and
  // removes go through `volume.*` and `particles.*`.
  const volumeNode = (v: VolumeRef, name: string): Node => {
    const ref = v.cell.transform.ref;
    let n = nodes.get(ref);
    if (!n) {
      n = new Node(ref, 'volume', 'volume', name, true, {
        visibility: (hide) => void api.run('volume.visibility', { volume: ref, visible: !hide }),
        remove: () => void api.run('volume.remove', { volume: ref }),
      });
      n.removable = true;
      nodes.set(ref, n);
    }
    const data = v.cell.obj?.data;
    const stats = data && volumeStats(data);
    n.set(name, stats ? `${gridText(stats)} · ${voxelText(stats)}` : undefined, !!v.cell.state.isHidden);
    return n;
  };

  const volumeReprNode = (r: VolumeRepresentationRef, volumeRef: string): Node => {
    const ref = r.cell.transform.ref;
    let n = nodes.get(ref);
    if (!n) {
      n = new Node(ref, 'volume-representation', 'volume', '', false, {
        visibility: (hide) =>
          void api.run('volume.visibility', { volume: volumeRef, representation: ref, visible: !hide }),
        remove: () => void api.run('volume.remove', { volume: volumeRef, representation: ref }),
      });
      n.removable = true;
      nodes.set(ref, n);
    }
    const iso = (r.cell.transform.params as { type?: { params?: { isoValue?: Volume.IsoValue } } } | undefined)?.type
      ?.params?.isoValue;
    n.set(
      r.cell.obj?.label ?? 'Representation',
      iso ? Volume.IsoValue.toString(iso) : undefined,
      !!r.cell.state.isHidden,
    );
    inspections.get(ref)?.form?.load(r.cell.transform.params as never);
    return n;
  };

  const particlesNode = (l: ParticleListRef, name: string): Node => {
    const ref = l.cell.transform.ref;
    let n = nodes.get(ref);
    if (!n) {
      n = new Node(ref, 'particles', 'shapes', name, true, {
        visibility: (hide) => void api.run('particles.visibility', { particles: ref, visible: !hide }),
        remove: () => void api.run('particles.remove', { particles: ref }),
      });
      n.removable = true;
      nodes.set(ref, n);
    }
    const count = l.cell.obj?.data.count;
    n.set(name, count !== undefined ? `${count.toLocaleString('en')} particles` : undefined, !!l.cell.state.isHidden);
    return n;
  };

  const particlesReprNode = (r: ParticleRepresentationRef, listRef: string): Node => {
    const ref = r.cell.transform.ref;
    let n = nodes.get(ref);
    if (!n) {
      n = new Node(ref, 'particles-representation', 'shapes', '', false, {
        visibility: (hide) =>
          void api.run('particles.visibility', { particles: listRef, representation: ref, visible: !hide }),
        remove: () => void api.run('particles.remove', { particles: listRef, representation: ref }),
      });
      n.removable = true;
      nodes.set(ref, n);
    }
    n.set(r.cell.obj?.label ?? 'Representation', undefined, !!r.cell.state.isHidden);
    inspections.get(ref)?.form?.load(r.cell.transform.params as never);
    return n;
  };

  const sync = () => {
    const list = plugin.managers.structure.hierarchy.current.structures;
    const names = structureNames(list);
    const seen = new Set<string>();
    const top: SceneNode[] = [];
    for (const s of list) {
      const ref = s.cell.transform.ref;
      const sn = structureNode(s, names.get(ref) ?? ref);
      seen.add(ref);
      const contents = contentsNode(s, ref);
      seen.add(contents.id);
      const kids: SceneNode[] = [contents];
      for (const c of s.components) {
        const cn = componentNode(c, ref);
        seen.add(cn.id);
        const reprs = c.representations.map((r) => {
          const rn = reprNode(r, cn.id, ref);
          seen.add(rn.id);
          return rn;
        });
        cn.setChildren(reprs);
        kids.push(cn);
      }
      sn.setChildren(kids);
      top.push(sn);
    }
    const volumes = plugin.managers.volume.hierarchy.current.volumes;
    const volumeLabels = namesOf(volumes);
    for (const v of volumes) {
      const ref = v.cell.transform.ref;
      const vn = volumeNode(v, volumeLabels.get(ref) ?? ref);
      seen.add(ref);
      vn.setChildren(
        v.representations.map((r) => {
          const rn = volumeReprNode(r, ref);
          seen.add(rn.id);
          return rn;
        }),
      );
      top.push(vn);
    }
    const lists = plugin.managers.particles.hierarchy.current.lists;
    const listLabels = namesOf(lists);
    for (const l of lists) {
      const ref = l.cell.transform.ref;
      const ln = particlesNode(l, listLabels.get(ref) ?? ref);
      seen.add(ref);
      ln.setChildren(
        l.representations.map((r) => {
          const rn = particlesReprNode(r, ref);
          seen.add(rn.id);
          return rn;
        }),
      );
      top.push(ln);
    }
    // Lazy contents' entity and chain nodes stay with their structure.
    for (const [id, n] of nodes) {
      const owner = id.includes('#') ? id.slice(0, id.indexOf('#')) : id;
      if (!seen.has(id) && !seen.has(owner)) {
        nodes.delete(id);
        inspections.delete(id);
        if (picked.value === n) picked.set(undefined);
      }
    }
    roots.set(top);
  };

  scope.add(plugin.managers.structure.hierarchy.behaviors.selection.subscribe(() => sync()));
  scope.add(plugin.managers.volume.hierarchy.behaviors.selection.subscribe(() => sync()));
  scope.add(plugin.managers.particles.hierarchy.behaviors.selection.subscribe(() => sync()));
  scope.add(
    plugin.state.data.events.cell.stateUpdated.subscribe((e) => {
      const n = nodes.get(e.ref);
      if (n) n.set(n.label.value, n.description.value, !!e.cell.state.isHidden);
    }),
  );

  return {
    roots: roots.asReadonly,
    picked: picked.asReadonly,
    pick: (n) => picked.set(n === picked.value ? undefined : n),
    inspect: (n) => {
      let inspection = inspections.get(n.id);
      if (!inspection) {
        inspection = inspect(plugin, api, n);
        inspections.set(n.id, inspection);
      }
      return inspection;
    },
    sync,
  };
}

function Structure_toLoci(s: Structure): StructureElement.Loci {
  return StructureElement.Loci.all(s);
}

function describeStructure(s: Structure): string {
  return `${s.elementCount.toLocaleString('en')} atoms`;
}

function reprLabel(plugin: PluginContext, r: StructureRepresentationRef): string {
  const name = (r.cell.transform.params as { type?: { name: string } } | undefined)?.type?.name;
  if (!name) return r.cell.obj?.label ?? 'Representation';
  return plugin.representation.structure.registry.get(name)?.label ?? name;
}

function colorLabel(plugin: PluginContext, r: StructureRepresentationRef): string | undefined {
  const name = (r.cell.transform.params as { colorTheme?: { name: string } } | undefined)?.colorTheme?.name;
  if (!name) return undefined;
  return plugin.representation.structure.themes.colorThemeRegistry.get(name)?.label ?? name;
}

/** Mol*'s component menu: add representation, color, modify by the selection, select this. */
function componentMenu(
  plugin: PluginContext,
  api: ViewerApi,
  ref: string,
  structureRef: string,
  loci: () => StructureElement.Loci | undefined,
): SceneMenuItem[] {
  const types = plugin.representation.structure.registry.types.map(([id, label]) => ({ id, label }));
  const themes = plugin.representation.structure.themes.colorThemeRegistry.types
    .slice(0, 24)
    .map(([id, label]) => ({ id, label }));
  return [
    {
      id: 'add-repr',
      label: 'Add representation',
      icon: 'plus',
      items: types.map((t) => ({
        id: t.id,
        label: t.label,
        run: () => void api.run('representation.set', { component: ref, structure: structureRef, type: t.id }),
      })),
    },
    {
      id: 'coloring',
      label: 'Set coloring',
      icon: 'palette',
      items: themes.map((t) => ({
        id: t.id,
        label: t.label,
        run: () => {
          const l = loci();
          if (l) void api.run('representation.color', { target: { loci: l }, theme: t.id });
        },
      })),
    },
    {
      id: 'modify',
      label: 'Modify by selection',
      icon: 'squares-unite',
      items: (['union', 'subtract', 'intersect'] as const).map((mode) => ({
        id: mode,
        label: mode === 'union' ? 'Include' : mode === 'subtract' ? 'Subtract' : 'Intersect',
        icon: mode === 'union' ? 'squares-unite' : mode === 'subtract' ? 'squares-subtract' : 'squares-intersect',
        run: () =>
          void api.run('component.modify', { component: ref, structure: structureRef, target: 'selection', mode }),
      })),
    },
    {
      id: 'select',
      label: 'Select this',
      icon: 'mouse-pointer-2',
      run: () => {
        const l = loci();
        if (l) void api.run('select.set', { target: { loci: l }, mode: 'set' });
      },
    },
    // Remove comes from the tree (the row is removable), set apart in the danger color.
  ];
}

/** A structure's entities, each with its chains (lazy contents, decision 18: a row's eye hides it). */
function entities(s: Structure, structureRef: string, nodes: Map<string, Node>, api: ViewerApi): SceneNode[] {
  const { entities } = s.model;
  const byEntity = new Map<string, Map<string, number>>();
  const loc = StructureElement.Location.create(s);
  for (const unit of s.units) {
    loc.unit = unit as Unit;
    loc.element = unit.elements[0];
    const entity = StructureProperties.chain.label_entity_id(loc);
    const chain = StructureProperties.chain.auth_asym_id(loc);
    let chains = byEntity.get(entity);
    if (!chains) byEntity.set(entity, (chains = new Map()));
    chains.set(chain, (chains.get(chain) ?? 0) + unit.elements.length);
  }
  const out: SceneNode[] = [];
  for (let i = 0; i < entities.data._rowCount; i++) {
    const id = entities.data.id.value(i);
    const chains = byEntity.get(id);
    if (!chains) continue;
    const key = `${structureRef}#entity:${id}`;
    let n = nodes.get(key);
    const label = entities.data.pdbx_description.value(i).join(', ') || `Entity ${id}`;
    if (!n) {
      n = new Node(key, 'entity', 'shapes', label, false, {
        visibility: (hide) =>
          void api.run('contents.visibility', { structure: structureRef, entity: id, visible: !hide }),
      });
      n.removable = false;
      nodes.set(key, n);
    }
    n.set(label, `${entities.data.type.value(i)} · ${chains.size} chain${chains.size === 1 ? '' : 's'}`, false);
    n.menu = contentsMenu(api, () => contentsLoci(s, id), label);
    const chainNodes: SceneNode[] = [];
    for (const [chain, atoms] of chains) {
      const ck = `${structureRef}#chain:${id}:${chain}`;
      let c = nodes.get(ck);
      if (!c) {
        c = new Node(ck, 'chain', 'sequence', `Chain ${chain}`, false, {
          visibility: (hide) => void api.run('contents.visibility', { structure: structureRef, chain, visible: !hide }),
        });
        c.removable = false;
        nodes.set(ck, c);
      }
      c.set(`Chain ${chain}`, `${atoms.toLocaleString('en')} atoms`, false);
      c.menu = contentsMenu(api, () => contentsLoci(s, id, chain), `${label} · chain ${chain}`);
      chainNodes.push(c);
    }
    n.setChildren(chainNodes);
    out.push(n);
  }
  return out;
}

/** The atoms of an entity, or of one of its chains, in a structure. */
function contentsLoci(s: Structure, entity: string, chain?: string): StructureElement.Loci {
  const isEntity = MS.core.rel.eq([MS.ammp('label_entity_id'), entity]);
  const test = chain ? MS.core.logic.and([isEntity, MS.core.rel.eq([MS.ammp('auth_asym_id'), chain])]) : isEntity;
  const query = MS.struct.generator.atomGroups({ 'chain-test': test });
  return StructureSelection.toLociWithSourceUnits(Script.getStructureSelection(query, s));
}

/** An entity's or a chain's row menu (parity.md §2: select, focus, a component of it). */
function contentsMenu(api: ViewerApi, loci: () => StructureElement.Loci, label: string): SceneMenuItem[] {
  return [
    {
      id: 'select',
      label: 'Select',
      icon: 'mouse-pointer-2',
      run: () => void api.run('select.set', { target: { loci: loci() }, mode: 'set' }),
    },
    { id: 'focus', label: 'Focus', icon: 'focus', run: () => void api.run('focus.set', { target: { loci: loci() } }) },
    {
      id: 'component',
      label: 'Make a component',
      icon: 'plus',
      run: () => void api.run('component.add', { target: { loci: loci() }, label }),
    },
  ];
}

/** What the inspector shows for a node: facts, and the parameters of a representation as a form. */
function inspect(plugin: PluginContext, api: ViewerApi, node: SceneNode): SceneInspection {
  const cell = plugin.state.data.cells.get(node.id);
  const slot = `inspector:${node.kind}`;
  if (node.kind === 'structure' && cell?.obj?.data) {
    const s = cell.obj.data as Structure;
    const model = s.model;
    const chains = new Set(s.units.map((u) => u.chainGroupId)).size;
    return {
      node,
      form: undefined,
      slot,
      facts: [
        ['Entry', model.entryId],
        ['Atoms', s.elementCount.toLocaleString('en')],
        ['Residues', s.polymerResidueCount.toLocaleString('en')],
        ['Chains', String(chains)],
        ['Model', model.label],
      ],
    };
  }
  if (node.kind === 'volume' && cell?.obj?.data) {
    const data = cell.obj.data as Volume;
    const stats = volumeStats(data);
    const signed = (x: number) => shortNumber(x).replace('-', '−');
    return {
      node,
      form: undefined,
      slot,
      facts: [
        ...(data.sourceData.kind ? [['Source', data.sourceData.kind.toUpperCase()] as const] : []),
        ['Grid', gridText(stats)],
        ['Voxel', voxelText(stats)],
        ['Mean · σ', `${signed(stats.mean)} · ${signed(stats.sigma)}`],
        ['Range', `${signed(stats.min)} … ${signed(stats.max)}`],
      ],
    };
  }
  if (node.kind === 'particles' && cell?.obj?.data) {
    const data = cell.obj.data as ParticleList;
    return {
      node,
      form: undefined,
      slot,
      facts: [
        ['Particles', data.count.toLocaleString('en')],
        ['Targets', String(data.targetInfo.size)],
      ],
    };
  }
  const updateOperation = UpdateOperations[node.kind];
  if (updateOperation && cell) {
    const parent = plugin.state.data.cells.get(cell.transform.parent);
    const def = cell.transform.transformer.definition.params?.(parent?.obj as never, plugin) as PD.Params | undefined;
    if (def) {
      let timer: ReturnType<typeof setTimeout> | undefined;
      const form = new ParamsFormModel(def, cell.transform.params as PD.Values<PD.Params>, {
        // Sliders change values at pointer rate: one update a moment after the last.
        onChange: (values) => {
          clearTimeout(timer);
          timer = setTimeout(() => void api.run(updateOperation, { representation: node.id, params: values }), 150);
        },
      });
      return { node, form, slot };
    }
  }
  return { node, form: undefined, slot };
}
