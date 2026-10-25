/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureElement } from '@molstar/model/model/structure';
import type { StructureComponentRef, StructureRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import type { Command, RepType } from '../language/ir.js';
import { RepTypes } from '../language/ir.js';
import { atomCount, compactExpression, componentLoci, polymerAtoms } from './atoms.js';
import { HostFailure, type Ops, type Where } from './ops.js';
import { MolstarReps, repTypeOf } from './reps.js';

// PyMOL's per-atom show flags on a tree of components (_spec §5.2). The tree is the truth: each command reads what
// every component of a structure shows, applies the command as set algebra on atoms, and writes back through the
// component operations. The console's own components (one per structure and type, in the state's `display`) grow and
// shrink; another component (a preset's, the person's) showing only hidden types is edited in place, and one showing
// other types too gives up the hidden type's representation, its remaining atoms moving to the console's component.

type Loci = StructureElement.Loci;

interface Shown {
  readonly component: StructureComponentRef;
  readonly ref: string;
  /** The console's own, and its type. */
  readonly own?: RepType;
  /** What it shows, in the IR's words (the console's own: its type); types PyMOL has no word for are left out. */
  readonly types: readonly RepType[];
  /** Mol* types it shows that PyMOL has no word for (carbohydrates, labels). */
  readonly other: readonly string[];
  readonly loci: Loci;
}

function shownOf(ops: Ops, s: StructureRef): Shown[] {
  const parent = s.cell.obj?.data;
  if (!parent) return [];
  const display = ops.holder.value.display;
  return s.components.map((component) => {
    const ref = component.cell.transform.ref;
    const own = display[ref];
    const types: RepType[] = [];
    const other: string[] = [];
    for (const r of component.representations) {
      const p = r.cell.transform.params as { type?: { name: string; params?: { visuals?: string[] } } } | undefined;
      const name = p?.type?.name ?? '';
      const t = own ?? repTypeOf(name, p?.type?.params);
      if (t) {
        if (!types.includes(t)) types.push(t);
      } else other.push(name);
    }
    return {
      component,
      ref,
      ...(own ? { own } : {}),
      types,
      other,
      loci: componentLoci(component.cell.obj?.data, parent),
    };
  });
}

const union = (xs: readonly Loci[], parent: StructureElement.Loci['structure']) =>
  xs.reduce((a, b) => StructureElement.Loci.union(a, b), StructureElement.Loci.none(parent));

/** The representations that draw polymers only. */
const PolymerTypes: ReadonlySet<RepType> = new Set(['cartoon', 'putty', 'tube', 'backbone']);

/** Adds atoms to the console's component of a type, making it when there is none. */
async function addOwn(
  ops: Ops,
  w: Where,
  shown: readonly Shown[],
  type: RepType,
  atoms: Loci,
): Promise<string | undefined> {
  // What another component shows as this type already: not drawn twice. A polymer representation takes the polymer
  // only (`show cartoon` after a load: the waters and the ligand left, PyMOL draws nothing more; Mol* would refuse).
  const others = shown.filter((x) => !x.own && x.types.includes(type)).map((x) => x.loci);
  const drawn = PolymerTypes.has(type) ? polymerAtoms(atoms) : atoms;
  const todo = others.length ? StructureElement.Loci.subtract(drawn, union(others, w.data)) : drawn;
  if (!atomCount(todo)) return undefined;
  const target = { molql: compactExpression(todo) };
  const own = shown.find((x) => x.own === type);
  if (own) {
    const r = await ops.run('component.modify', { component: own.ref, target, mode: 'add', structure: w.ref });
    return r.did;
  }
  const rep = MolstarReps[type];
  const r = await ops.run(
    'component.add',
    { target, label: type, representation: rep.type, structure: w.ref },
    { copy: true },
  );
  const created = (r.value as readonly { ref: string }[] | undefined)?.[0];
  if (!created) throw new HostFailure(`Could not make a component for ${type}.`);
  ops.holder.update({ display: { ...ops.holder.value.display, [created.ref]: type } });
  await styleOwn(ops, created.ref, type);
  return r.did;
}

/** A new component's representation, styled as the type asks (licorice sticks, a wireframe mesh, a tube). */
async function styleOwn(ops: Ops, componentRef: string, type: RepType): Promise<void> {
  const rep = MolstarReps[type];
  if (!rep.typeParams && !rep.sizeTheme) return;
  const s = ops.structures().find((x) => x.components.some((c) => c.cell.transform.ref === componentRef));
  const c = s?.components.find((x) => x.cell.transform.ref === componentRef);
  const r = c?.representations[0];
  if (!r) return;
  await ops.run('representation.update', {
    node: r.cell.transform.ref,
    params: {
      ...(rep.typeParams ? { type: { name: rep.type, params: rep.typeParams } } : {}),
      ...(rep.sizeTheme ? { sizeTheme: rep.sizeTheme } : {}),
    },
  });
}

/** Takes atoms out of what shows a type: the console's component, and the others (in place, or by splitting). */
async function hideType(ops: Ops, w: Where, type: RepType, atoms: Loci): Promise<number> {
  let changed = 0;
  for (const x of shownOf(ops, w.structure)) {
    if (!x.types.includes(type) || !StructureElement.Loci.areIntersecting(x.loci, atoms)) continue;
    const target = { molql: compactExpression(StructureElement.Loci.intersect(x.loci, atoms)) };
    if (x.own || (x.types.length === 1 && !x.other.length)) {
      // Shows only this type: the atoms leave the component (Prism removes one left empty).
      await ops.run('component.modify', { component: x.ref, target, mode: 'remove', structure: w.ref });
    } else {
      // Shows other types too: it gives up this one; its atoms not hidden keep it in the console's component.
      const keep = StructureElement.Loci.subtract(x.loci, atoms);
      for (const r of x.component.representations) {
        const p = r.cell.transform.params as { type?: { name: string; params?: { visuals?: string[] } } } | undefined;
        if (repTypeOf(p?.type?.name ?? '', p?.type?.params) === type) {
          await ops.run('representation.remove', { node: r.cell.transform.ref });
        }
      }
      if (atomCount(keep)) await addOwn(ops, w, shownOf(ops, w.structure), type, keep);
    }
    changed++;
  }
  return changed;
}

/** The types a structure shows now (for `hide everything`). */
function typesShown(shown: readonly Shown[]): RepType[] {
  return RepTypes.filter((t) => shown.some((x) => x.types.includes(t)));
}

/** `show`, `hide`, `as`: per structure the selection selects atoms in. */
export async function display(ops: Ops, cmd: Extract<Command, { kind: 'display' }>): Promise<string> {
  const where = cmd.mode === 'hide' ? ops.where(cmd.sel) : ops.some(cmd.sel);
  let atoms = 0;
  for (const w of where) {
    const shown = shownOf(ops, w.structure);
    const hide = cmd.mode === 'show' ? [] : cmd.mode === 'only' || cmd.reps === 'all' ? typesShown(shown) : cmd.reps;
    for (const t of hide) {
      if (cmd.mode === 'only' && cmd.reps !== 'all' && cmd.reps.includes(t)) continue;
      await hideType(ops, w, t, w.loci);
    }
    // Everything: also what PyMOL has no word for (carbohydrate symbols, labels), which leaves with the atoms.
    if (cmd.mode === 'only' || (cmd.mode === 'hide' && cmd.reps === 'all')) {
      for (const x of shownOf(ops, w.structure)) {
        if (!x.other.length || !StructureElement.Loci.areIntersecting(x.loci, w.loci)) continue;
        const target = { molql: compactExpression(StructureElement.Loci.intersect(x.loci, w.loci)) };
        await ops.run('component.modify', { component: x.ref, target, mode: 'remove', structure: w.ref });
      }
    }
    if (cmd.mode !== 'hide') {
      if (cmd.reps === 'all')
        throw new HostFailure('Showing everything at once is not supported.', 'Name a representation.');
      for (const t of cmd.reps) await addOwn(ops, w, shownOf(ops, w.structure), t, w.loci);
    }
    atoms += atomCount(w.loci);
  }
  const what = cmd.reps === 'all' ? 'everything' : cmd.reps.join(', ');
  const verb = cmd.mode === 'hide' ? 'Hid' : 'Showed';
  return `${verb} ${what} for ${atoms} atom${atoms === 1 ? '' : 's'}${where.length > 1 ? ` in ${where.length} structures` : ''}.`;
}

/** `cartoon putty|tube|automatic`: what cartoon-like atoms of the selection show as. */
export async function cartoon(ops: Ops, cmd: Extract<Command, { kind: 'cartoon' }>): Promise<string> {
  const from: readonly RepType[] =
    cmd.style === 'automatic'
      ? ['putty', 'tube']
      : (['cartoon', 'putty', 'tube'].filter((t) => t !== cmd.style) as RepType[]);
  const to: RepType = cmd.style === 'automatic' ? 'cartoon' : cmd.style;
  let atoms = 0;
  for (const w of ops.some(cmd.sel)) {
    const shown = shownOf(ops, w.structure);
    const drawn = union(
      shown.filter((x) => from.some((t) => x.types.includes(t))).map((x) => x.loci),
      w.data,
    );
    const y = StructureElement.Loci.intersect(drawn, w.loci);
    if (!atomCount(y)) continue;
    for (const t of from) await hideType(ops, w, t, y);
    await addOwn(ops, w, shownOf(ops, w.structure), to, y);
    atoms += atomCount(y);
  }
  if (!atoms) return `Nothing in '${cmd.sel.text}' is shown as a cartoon: show it first (show cartoon, …).`;
  return `Drew the cartoon of ${atoms} atoms as ${to}.`;
}
