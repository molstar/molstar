/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A section's blocks (_spec/README.md "Resolved view", Blocks): its `block` children, else a complex's units (a section
// `of` a complex: ordered by `classes`, polymers in the complex's order, then each class, EC-101); each block's form
// (`whole`, `residues`, `atoms`), its label (a chip's: its molecule's name, else its one component, else its id), its
// direction and scale, its copies merged (EC-65), its group (an assembly of several); the section's axis over them.

import { exactAddress, type Node, printAddress, type Value } from '@molstar/core/repr/tree';
import {
  type AssemblyEntry,
  atomsSystemOf,
  describedAddress,
  type StaticModel,
  type UnitEntry,
} from '@molstar/seqx-core/resolve';
import { axisTotal } from './axis.js';
import type { Axis, BlockForm, ResolvedBlock } from './view.js';

const param = (n: Node, key: string): Value | undefined => n.params?.[key];

/** The classes a section of a complex shows when it names none: water left out (EC-101). */
export const DefaultSectionClasses: readonly string[] = ['polymer', 'branched', 'ligand', 'ion'];

/** A section's `of` as an address: an exact one, or a node a spec source holds (`<source>/(node :key <key>)`). */
function sectionOf(model: StaticModel, section: Node): string | undefined {
  return refOf(model, param(section, 'of'));
}

/** The complex a section is `of`, when it is one (a source's while held: an opened model's). */
export function sectionComplex(model: StaticModel, section: Node): string | undefined {
  const of = sectionOf(model, section);
  return of !== undefined && model.molecular.complexes.has(of) ? of : undefined;
}

/**
 * The source a section's `of` names a complex of when that complex is a node the source holds, or will
 * (`<source>/(node :key <key>)`, the source the spec's): read while it is held (an opened model), pending before
 * (_spec/README.md "Resolved view", A section of a source's complex).
 */
export function openedComplex(model: StaticModel, section: Node): string | undefined {
  const of = sectionOf(model, section);
  const at = of?.indexOf('/(node ') ?? -1;
  const source = at > 0 ? of!.slice(0, at) : undefined;
  return source !== undefined && model.sources.has(source) ? source : undefined;
}

const Opens = new WeakMap<StaticModel, boolean>();

/**
 * Whether a view of the spec has a section `of` a source's complex (`openedComplex`): what resolves its views over the
 * opened model (`openedModel`); kept per model. A spec without one never reads the overlay.
 */
export function opensComplexes(model: StaticModel): boolean {
  let opens = Opens.get(model);
  if (opens !== undefined) return opens;
  opens = false;
  for (const c of model.spec.children ?? [])
    for (const v of c.kind === 'views' ? (c.children ?? []) : [])
      for (const s of v.children ?? []) if (s.kind === 'section' && openedComplex(model, s)) opens = true;
  Opens.set(model, opens);
  return opens;
}

/** The assemblies a section is `of` (one, or a list: a group each), when it is of assemblies. */
export function sectionAssemblies(model: StaticModel, section: Node): AssemblyEntry[] | undefined {
  const of = param(section, 'of');
  const list = (Array.isArray(of) ? of : [of]).map((v) => model.assemblies.get(exactAddress(v) ?? ''));
  return list.length && list.every((a) => a !== undefined) ? (list as AssemblyEntry[]) : undefined;
}

/**
 * The system a block's `of` names: an exact address, or a derived sequence's (`/data/hg38/(sequence :name chr1)`, a
 * form no exact address reads).
 */
export function systemRef(model: StaticModel, of: Value | undefined): string | undefined {
  return typeof of === 'string' && model.systems.has(of) ? of : refOf(model, of);
}

/** An address as written: an exact one, or a node a spec source holds (`<source>/(node :key <key>)…`), canonical. */
function refOf(model: StaticModel, v: Value | undefined): string | undefined {
  return exactAddress(v) ?? describedAddress(v, (a) => model.sources.has(a));
}

/**
 * The overrides of a section of a source's complex that name no unit an opened model over this spec has held (a
 * mistyped id): their ids, to be said; an override of a unit an update dropped is not one.
 */
export function unknownOverrides(model: StaticModel, section: Node): string[] {
  const complex = sectionComplex(model, section);
  if (complex === undefined || openedComplex(model, section) === undefined) return [];
  const held = model.heldUnits;
  const out: string[] = [];
  for (const c of section.children ?? []) {
    if (c.kind !== 'block' || c.id === undefined || !isOverride(model, c, complex)) continue;
    const address = `${complex}/${c.id}`;
    if (!model.molecular.units.has(address) && !held?.has(address)) out.push(c.id);
  }
  return out;
}

/** A `block` with no `of`, or one naming its id's unit of the complex: an override. */
const isOverride = (model: StaticModel, n: Node, complex: string | undefined) => {
  const of = param(n, 'of');
  return of === undefined || refOf(model, of) === `${complex}/${n.id}`;
};

/** A derived block's address: the keyed form under its section, never an id step (it has no id). */
const DerivedAddress = new WeakMap<Node, string>();
/** A derived block's sequence's printed name (the `system` field's value). */
const DerivedName = new WeakMap<Node, string>();
/** A derived block's group: its assembly's index in the section's `of` (one assembly one group; a list a group each, EC-210). */
const DerivedGroup = new WeakMap<Node, number>();

/** A group of a section of assemblies' blocks: its assembly, its label (the assembly's name), its blocks' run. */
export interface BlockGroup {
  readonly assembly: string;
  readonly label: string;
  /** Its first block's index in the section's blocks, and how many it has. */
  readonly first: number;
  readonly count: number;
}

/** The derived systems a section of assemblies' `block` child overrides: its `of`'s system, else each its id names. */
function overrideSystems(model: StaticModel, c: Node, assemblies: readonly AssemblyEntry[]): string[] {
  const of = systemRef(model, param(c, 'of'));
  if (of !== undefined) return [of];
  const out: string[] = [];
  if (c.id !== undefined)
    for (const a of assemblies) {
      const s = a.systemAt[a.find(c.id) ?? -1];
      if (s) out.push(s);
    }
  return out;
}

/**
 * The `block` children of a section of assemblies that override a derived system (_spec/README.md "Resolved view",
 * Blocks): one whose `of` names it (an exact address or a sequence's form, an alias too), or with no `of` whose id names
 * it (a name or an alias by the assembly's names rule); the rule `assemblyBlocks` applies, so an editor finds the
 * override already written, whatever its spelling. Empty when the section is not of assemblies.
 */
export function blockOverridesOf(model: StaticModel, section: Node, system: string): Node[] {
  const assemblies = sectionAssemblies(model, section);
  if (!assemblies) return [];
  return (section.children ?? []).filter(
    (c) => c.kind === 'block' && overrideSystems(model, c, assemblies).includes(system),
  );
}

/**
 * A section of assemblies' blocks (_spec/README.md "Resolved view", Blocks): a block a derived system, each assembly's
 * in its order and after its `include`, then the section's `include`; a `block` child whose `of` is a derived system
 * (or whose id is its name) overrides that block only (`direction`, `scale`, `label`, `range`). At most `max` blocks;
 * `total` is how many there were.
 */
export function assemblyBlocks(
  model: StaticModel,
  section: Node,
  assemblies: readonly AssemblyEntry[],
  at: string,
): { nodes: Node[]; total: number } {
  const over = new Map<string, Record<string, Value>>();
  for (const c of section.children ?? []) {
    if (c.kind !== 'block') continue;
    const params: Record<string, Value> = {};
    for (const k of ['direction', 'scale', 'label', 'range', 'form'] as const)
      if (param(c, k) !== undefined) params[k] = param(c, k)!;
    for (const s of overrideSystems(model, c, assemblies)) over.set(s, params);
  }
  const nodes: Node[] = [];
  let total = 0;
  const max = model.config.maxAssemblyBlocks;
  for (let g = 0; g < assemblies.length; g++) {
    const a = assemblies[g];
    const keep = a.include(param(section, 'include'));
    for (const i of a.order) {
      const system = a.systemAt[i];
      const s = a.sequences[i];
      if (!system || (keep && !keep(i, s))) continue;
      total++;
      if (nodes.length >= max) continue;
      const label = model.systems.get(system)?.label ?? s.name;
      const n: Node = { kind: 'block', params: { of: system, label, ...over.get(system) } };
      DerivedAddress.set(
        n,
        printAddress([...at.split('/').filter(Boolean), { kind: 'block', params: { name: s.name } }]),
      );
      DerivedName.set(n, s.name);
      DerivedGroup.set(n, g);
      nodes.push(n);
    }
  }
  return { nodes, total };
}

/**
 * A section's block nodes: its `block` children (explicit blocks win), else a block a unit of the complex it is `of`,
 * by `classes` (a class not listed left out), polymers in the complex's order, then each class in order. A `block`
 * child whose id is a derived unit's (no `of`, or that unit) overrides that unit only, its `form` and `label` (a
 * merged copies block by its first unit's id: every copy); a section whose children are all overrides lists every
 * derived unit, as units are added too.
 */
export function blockNodes(model: StaticModel, section: Node, at = ''): Node[] {
  const assemblies = sectionAssemblies(model, section);
  if (assemblies) return assemblyBlocks(model, section, assemblies, at).nodes;
  const own = (section.children ?? []).filter((c) => c.kind === 'block' && c.id !== undefined);
  const complex = sectionComplex(model, section);
  const units = complex === undefined ? undefined : model.molecular.complexes.get(complex)?.units;
  // a source's complex not held (yet, or any more): no block, its overrides waiting
  const opened = openedComplex(model, section) !== undefined;
  if (!units) return opened ? [] : own;
  const byId = new Map(units.map((u) => [u.id, u]));
  // over a source's complex a block naming no other system overrides its unit, held or not: a unit an update drops
  // takes its override with it, never turning the rest into a list
  const override = (n: Node) => (byId.has(n.id!) || opened) && isOverride(model, n, complex);
  if (own.length && !own.every(override)) return own;
  // each override's params by the units it reaches: its unit, or with merged copies its unit's copy group
  const merge = param(section, 'copies') !== 'each';
  const over = new Map<string, Record<string, Value>>();
  for (const n of own) {
    const u = byId.get(n.id!);
    if (!u) continue;
    const params: Record<string, Value> = {};
    for (const k of ['form', 'label'] as const) if (param(n, k) !== undefined) params[k] = param(n, k)!;
    for (const v of units) if (v === u || (merge && v.copy === u.copy)) over.set(v.id, params);
  }
  const written = param(section, 'classes');
  const classes = Array.isArray(written) ? written.map(String) : DefaultSectionClasses;
  const out: Node[] = [];
  for (const c of classes)
    for (const u of units)
      if (u.class === c) out.push({ kind: 'block', id: u.id, params: { of: u.address, ...over.get(u.id) } });
  return out;
}

/**
 * A block's display units in its section: one for a `whole` block (a chip: its molecule one unit, _spec/README.md
 * "Render", Section layouts), else its positions shown; undefined while its length is unknown.
 */
export function blockUnits(block: { readonly form?: string; readonly axis: Axis }): number | undefined {
  if (block.form === 'whole' && block.axis.ranges) return 1;
  return axisTotal(block.axis);
}

/** A section's blocks before their layers, and the section's axis over them. */
export interface BlockAxes {
  readonly axis: Axis;
  readonly blocks: readonly Omit<ResolvedBlock, 'layers'>[];
  /** A section of several assemblies' groups, in order (EC-210). */
  readonly groups?: readonly BlockGroup[];
}

/** A block's scale as written: a positive number, else 1 (_spec/README.md "Resolved view", Blocks). */
export function blockScale(node: Node): number {
  const s = param(node, 'scale');
  return typeof s === 'number' && s > 0 && Number.isFinite(s) ? s : 1;
}

type AxisOf = (node: Node, at: string, system: string, dir: 1 | -1) => Axis;
type Report = (message: string, at: string) => void;

/**
 * A block's form and the system it shows (_spec/README.md "Resolved view", Blocks): `form` as written, else `whole`
 * for a molecule or a unit of one and `residues` for the rest; `atoms` shows the molecule's (or unit's) atoms system,
 * its range the atoms of the residues in range; refused (the default form) where there is no atoms system.
 */
function formOf(model: StaticModel, node: Node, system: string): { form: BlockForm; system: string; node: Node } {
  const m = model.molecular;
  const unit = m.units.get(system);
  const molecule = m.molecules.get(unit ? unit.of : system);
  const written = param(node, 'form');
  const fallback: BlockForm = molecule ? 'whole' : 'residues';
  if (written === 'whole' || written === 'residues') return { form: written, system, node };
  if (written !== 'atoms') return { form: fallback, system, node };
  const atoms = m.atoms.get(atomsSystemOf(system));
  if (!atoms) return { form: fallback, system, node };
  const range = param(node, 'range');
  if (!Array.isArray(range) || range.length !== 2) return { form: 'atoms', system: atoms.address, node };
  // the atoms of the residues in range: a residue's atoms are a run (`AtomsEntry.offsets`)
  const o = atoms.offsets;
  const n = o.length - 1;
  const clamp = (r: number) => o[Math.max(0, Math.min(n, Number(r)))];
  const params = { ...node.params, range: [clamp(range[0] as number), clamp(range[1] as number)] };
  return { form: 'atoms', system: atoms.address, node: { ...node, params } };
}

/** A block's label: as written, else its unit's, else (a molecule) its name or its one component, else its id. */
function labelOf(model: StaticModel, node: Node, system: string, unit: UnitEntry | undefined): string {
  const written = param(node, 'label');
  if (typeof written === 'string') return written;
  if (unit?.label !== undefined) return unit.label;
  const molecule = model.molecular.molecules.get(unit ? unit.of : system);
  if (molecule) {
    if (molecule.name !== undefined) return molecule.name;
    const r = molecule.residues;
    const start = unit?.start ?? 0;
    const end = unit?.end ?? r?.length ?? 0;
    if (r && end - start === 1 && r[start]) return r[start];
  }
  return node.id!;
}

/**
 * A section's blocks (_spec/README.md "Resolved view", Blocks): each its axis (its system, its `range`, the section's
 * direction times its own `direction`) and form, its `scale` (its display units its positions times it); copies (units
 * of one copy group, or the same system, with the same range, form, direction and scale) merged into the first's place
 * with `copies: 'merge'` (the default), labeled `×N`; the section's axis their display units in order: on one system in
 * a row (or on a ring), each block the section's direction at scale 1, the blocks' pieces end to end, else units on no
 * system (the section's address, `pieces: 'blocks'`), their total the sum in a row and the longest block in a stack;
 * a section of several assemblies' groups. A block with no system is reported and left out.
 */
export function blockAxes(
  model: StaticModel,
  section: Node,
  nodes: readonly Node[],
  at: string,
  dir: 1 | -1,
  ring: boolean,
  axisOf: AxisOf,
  report: Report,
): BlockAxes {
  const merge = param(section, 'copies') !== 'each';
  const assemblies = sectionAssemblies(model, section);
  const assemblyOf = (g: number) => assemblies?.[g]?.address ?? '';
  const assemblyLabel = (g: number) => assemblies?.[g]?.names[0] ?? assemblyOf(g);
  type Group = {
    node: Node;
    address: string;
    axis: Axis;
    key: string;
    members: string[];
    systems: string[];
    form: BlockForm;
    label: string;
    dir: 1 | -1;
    scale: number;
    group: number | undefined;
    name: string;
  };
  const groups: Group[] = [];
  // copy groups by key: linear in the blocks (an assembly's thousands of contigs)
  const byKey = new Map<string, Group>();
  for (const n of nodes) {
    const address = DerivedAddress.get(n) ?? `${at}/${n.id}`;
    const of = systemRef(model, param(n, 'of'));
    if (of === undefined || !model.systems.has(of)) {
      report('The block shows no system of the spec', address);
      continue;
    }
    const unit = model.molecular.units.get(of);
    const { form, system, node } = formOf(model, n, of);
    const range = param(node, 'range');
    // a block's direction against the section's, its scale (_spec/README.md "Resolved view", Blocks)
    const own: 1 | -1 = param(n, 'direction') === 'reversed' ? (-dir as 1 | -1) : dir;
    const scale = blockScale(n);
    // units of one copy group are copies whatever their systems (EC-65): equal `of`, range, direction and scale
    const key = `${unit ? unit.copy : system}|${form}|${Array.isArray(range) ? range.join(',') : 'all'}${own === dir ? '' : '|r'}${scale === 1 ? '' : `|x${scale}`}`;
    const copy = merge ? byKey.get(key) : undefined;
    if (copy) {
      copy.members.push(address);
      copy.systems.push(system);
      continue;
    }
    const axis = axisOf(node, address, system, own);
    const label = labelOf(model, n, of, unit);
    const group = DerivedGroup.get(n);
    const name = DerivedName.get(n) ?? model.systems.get(of)?.name ?? of;
    const g = {
      node,
      address,
      axis,
      key,
      members: [address],
      systems: [system],
      form,
      label,
      dir: own,
      scale,
      group,
      name,
    };
    groups.push(g);
    if (!byKey.has(key)) byKey.set(key, g);
  }
  const stack = !ring && param(section, 'arrange') === 'stack';
  let offset: number | undefined = 0;
  let longest = 0;
  // a section of assemblies: a group each with blocks, its blocks a run in order (EC-210); a block's `group` its index
  // in that list (an assembly that `include` leaves without blocks has no group)
  let groupList: BlockGroup[] | undefined;
  const groupOf: (number | undefined)[] = [];
  for (let k = 0; k < groups.length; k++) {
    const gi = groups[k].group;
    if (gi === undefined) continue;
    groupList ??= [];
    const last = groupList[groupList.length - 1];
    if (last && last.assembly === assemblyOf(gi)) (last as { count: number }).count++;
    else groupList.push({ assembly: assemblyOf(gi), label: assemblyLabel(gi), first: k, count: 1 });
    groupOf[k] = groupList.length - 1;
  }
  const blocks = groups.map((g, k) => {
    // a chip is one unit, however many residues it stands for
    const units = blockUnits(g);
    // a chip (a `whole` block, one position) is one unit at any scale: its slot is fixed
    if (g.form === 'whole' || units === 1) g.scale = 1;
    const total = units === undefined ? undefined : units * g.scale;
    const start = stack ? 0 : offset;
    if (offset !== undefined) offset = total === undefined ? undefined : offset + total;
    if (total !== undefined) longest = Math.max(longest, total);
    return {
      address: g.address,
      label: g.members.length > 1 ? `${g.label} ×${g.members.length}` : g.label,
      axis: g.axis,
      form: g.form,
      ...(start !== undefined ? { offset: start } : {}),
      members: g.members,
      systems: g.systems,
      copy: g.key,
      name: g.name,
      ...(g.scale !== 1 ? { scale: g.scale } : {}),
      ...(groupOf[k] !== undefined ? { group: groupOf[k] } : {}),
    };
  });
  const first = groups[0]?.axis;
  const system = first?.system ?? '';
  const pending = groups.some((g) => !g.axis.ranges);
  // one system's blocks are one axis only when they read one way, the section's, at their own length (genome-coordinates:
  // blocks on one system merge only on one direction)
  if (
    first &&
    !stack &&
    groups.every((g) => g.axis.system === system && g.form !== 'whole' && g.dir === dir && g.scale === 1)
  ) {
    // pieces of one sequence: one axis, its breaks between the blocks
    const ranges = pending ? undefined : groups.flatMap((g) => g.axis.ranges!);
    return {
      axis: {
        system,
        ...(first.length !== undefined ? { length: first.length } : {}),
        cyclic: first.cyclic,
        dir,
        pieces: 'blocks',
        ...(ranges ? { ranges } : {}),
      },
      blocks,
      ...(groupList ? { groups: groupList } : {}),
    };
  }
  // units of the section, on no system: a position on them is no block's (a selection, a focus, a bookmark or a cursor
  // on a system never lands on the wrong block's); a block's hits carry its own system
  const total = stack ? longest : offset;
  return {
    axis: {
      system: at,
      ...(pending || total === undefined ? {} : { length: total, ranges: [0, total] }),
      cyclic: false,
      dir: 1,
      pieces: 'blocks',
    },
    blocks,
    ...(groupList ? { groups: groupList } : {}),
  };
}
