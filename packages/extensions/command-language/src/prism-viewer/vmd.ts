/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color } from '@molstar/core/util/color';
import { StructureElement } from '@molstar/model/model/structure';
import type {
  StructureComponentRef,
  StructureRef,
  StructureRepresentationRef,
} from '@molstar/plugin/state/manager/structure/hierarchy-state';
import type { Command, ColorSpec, VmdRep } from '../language/ir.js';
import { selectionOn, vmdSelection } from '../language/selection.js';
import { atomCount, compactExpression, componentLoci, lociOf } from './atoms.js';
import { HostFailure, type Ops } from './ops.js';
import { colorTheme, MolstarReps } from './reps.js';

// VMD's reps on the tree (_spec §5.3): a molecule's reps are its components' representations in tree order, so a
// structure loaded any way has reps, and Undo keeps them in step with the tree. A rep the console adds is a component
// of its selection with one representation, whose own params carry the rep's coloring and material.

interface RepEntry {
  readonly structure: StructureRef;
  readonly component: StructureComponentRef;
  readonly representation: StructureRepresentationRef;
}

function repsOf(s: StructureRef): RepEntry[] {
  return s.components.flatMap((component) =>
    component.representations.map((representation) => ({ structure: s, component, representation })),
  );
}

function entry(ops: Ops, object: string, index: number): RepEntry {
  const s = ops.structure(ops.ref(object));
  if (!s) throw new HostFailure(`There is no molecule '${object}'.`);
  const reps = repsOf(s);
  const e = reps[index];
  if (!e) {
    throw new HostFailure(
      `${object} has no rep ${index}.`,
      reps.length ? `Its reps are 0 to ${reps.length - 1}.` : 'It has no reps: add one with mol addrep.',
    );
  }
  return e;
}

const TransparentAlpha = 0.5;

/** The params of a rep's representation that say its coloring and material. */
function repParams(
  rep: Pick<VmdRep, 'color' | 'material'>,
  type: string,
  typeParams?: Readonly<Record<string, unknown>>,
) {
  const theme: { name: string; params: Record<string, unknown> } =
    'rgb' in rep.color ? { name: 'uniform', params: { value: Color(rep.color.rgb) } } : colorTheme(rep.color);
  return {
    colorTheme: theme,
    type: {
      name: type,
      params: { ...(typeParams ?? {}), alpha: rep.material === 'transparent' ? TransparentAlpha : 1 },
    },
  };
}

export async function addRep(ops: Ops, cmd: Extract<Command, { kind: 'rep.add' }>): Promise<string> {
  const ref = ops.ref(cmd.object);
  const s = ops.structure(ref);
  const data = s?.cell.obj?.data;
  if (!s || !data) throw new HostFailure(`${cmd.object} is not built yet.`);
  const sel = vmdSelection(cmd.rep.selection, cmd.object, cmd.at);
  if (!atomCount(lociOf(sel, cmd.object, data))) {
    throw new HostFailure(
      `'${cmd.rep.selection}' selects nothing in ${cmd.object}.`,
      'VMD would add an empty rep; Mol* adds none. Check the selection (mol selection …).',
    );
  }
  const rep = MolstarReps[cmd.rep.style];
  const before = new Set(s.components.map((c) => c.cell.transform.ref));
  await ops.run(
    'component.add',
    {
      target: { molql: selectionOn(sel, cmd.object) },
      label: `${cmd.rep.style}: ${cmd.rep.selection}`,
      representation: rep.type,
      structure: ref,
    },
    { copy: true },
  );
  const made = ops.structure(ref)?.components.find((c) => !before.has(c.cell.transform.ref));
  const r = made?.representations[0];
  if (!r) throw new HostFailure(`Could not add the rep to ${cmd.object}.`);
  await ops.run('representation.update', {
    node: r.cell.transform.ref,
    params: { ...repParams(cmd.rep, rep.type, rep.typeParams), ...(rep.sizeTheme ? { sizeTheme: rep.sizeTheme } : {}) },
  });
  return `Added rep ${repsOf(ops.structure(ref)!).length - 1} to ${cmd.object}: ${cmd.rep.style} of '${cmd.rep.selection}'.`;
}

export async function removeRep(ops: Ops, cmd: Extract<Command, { kind: 'rep.remove' }>): Promise<string> {
  const e = entry(ops, cmd.object, cmd.index);
  if (e.component.representations.length === 1) {
    await ops.run('component.remove', {
      component: e.component.cell.transform.ref,
      structure: e.structure.cell.transform.ref,
    });
  } else {
    await ops.run('representation.remove', { node: e.representation.cell.transform.ref });
  }
  return `Removed rep ${cmd.index} of ${cmd.object}.`;
}

export async function showRep(ops: Ops, cmd: Extract<Command, { kind: 'rep.visibility' }>): Promise<string> {
  const e = entry(ops, cmd.object, cmd.index);
  await ops.run('representation.visibility', { node: e.representation.cell.transform.ref, visible: cmd.visible });
  return `${cmd.visible ? 'Showed' : 'Hid'} rep ${cmd.index} of ${cmd.object}.`;
}

type RepTransformParams = {
  type?: { name: string; params?: Record<string, unknown> };
  colorTheme?: { name: string; params?: Record<string, unknown> };
};

export async function modifyRep(ops: Ops, cmd: Extract<Command, { kind: 'rep.modify' }>): Promise<string> {
  const e = entry(ops, cmd.object, cmd.index);
  // A rep that shares its component with others (a preset's): its selection is theirs, and a new style would move it
  // to the component's end (the reps after it would change their indices).
  if ((cmd.patch.selection !== undefined || cmd.patch.style !== undefined) && e.component.representations.length > 1) {
    throw new HostFailure(
      `Rep ${cmd.index} of ${cmd.object} shares its selection with other reps: changing it would change them too.`,
      'Add a rep of its own instead: mol addrep, then mol delrep this one.',
    );
  }
  const ref = e.structure.cell.transform.ref;
  const node = e.representation.cell.transform.ref;
  const current = e.representation.cell.transform.params as RepTransformParams | undefined;
  const changed: string[] = [];
  if (cmd.patch.selection !== undefined) {
    const data = e.structure.cell.obj?.data;
    if (!data) throw new HostFailure(`${cmd.object} is not built yet.`);
    const sel = vmdSelection(cmd.patch.selection, cmd.object, cmd.at);
    const want = lociOf(sel, cmd.object, data);
    if (!atomCount(want)) throw new HostFailure(`'${cmd.patch.selection}' selects nothing in ${cmd.object}.`);
    const has = componentLoci(e.component.cell.obj?.data, data);
    const add = StructureElement.Loci.subtract(want, has);
    const remove = StructureElement.Loci.subtract(has, want);
    const component = e.component.cell.transform.ref;
    if (atomCount(add))
      await ops.run('component.modify', {
        component,
        target: { molql: compactExpression(add) },
        mode: 'add',
        structure: ref,
      });
    if (atomCount(remove))
      await ops.run('component.modify', {
        component,
        target: { molql: compactExpression(remove) },
        mode: 'remove',
        structure: ref,
      });
    changed.push('selection');
  }
  if (cmd.patch.style !== undefined) {
    const rep = MolstarReps[cmd.patch.style];
    const alpha = current?.type?.params?.alpha;
    const component = e.component.cell.transform.ref;
    // A new type: the representation is made again with the rep's coloring.
    await ops.run('representation.remove', { node });
    await ops.run('representation.set', { component, type: rep.type, params: rep.typeParams ?? {}, structure: ref });
    const made = ops
      .structure(ref)
      ?.components.find((c) => c.cell.transform.ref === component)
      ?.representations.find(
        (r) => (r.cell.transform.params as RepTransformParams | undefined)?.type?.name === rep.type,
      );
    if (made) {
      await ops.run('representation.update', {
        node: made.cell.transform.ref,
        params: {
          ...(current?.colorTheme ? { colorTheme: current.colorTheme } : {}),
          ...(alpha !== undefined ? { type: { name: rep.type, params: { ...(rep.typeParams ?? {}), alpha } } } : {}),
          ...(rep.sizeTheme ? { sizeTheme: rep.sizeTheme } : {}),
        },
      });
    }
    changed.push('style');
    return done(cmd, changed);
  }
  if (cmd.patch.color !== undefined || cmd.patch.material !== undefined) {
    const color: ColorSpec | undefined = cmd.patch.color;
    const type = current?.type?.name ?? '';
    const params: Record<string, unknown> = {};
    if (color)
      params.colorTheme = 'rgb' in color ? { name: 'uniform', params: { value: Color(color.rgb) } } : colorTheme(color);
    if (cmd.patch.material !== undefined) {
      // The type's params as they are, the alpha changed (an update replaces the type's params whole).
      params.type = {
        name: type,
        params: {
          ...(current?.type?.params ?? {}),
          alpha: cmd.patch.material === 'transparent' ? TransparentAlpha : 1,
        },
      };
    }
    await ops.run('representation.update', { node, params });
    if (color) changed.push('coloring');
    if (cmd.patch.material !== undefined) changed.push('material');
  }
  return done(cmd, changed);
}

const done = (cmd: Extract<Command, { kind: 'rep.modify' }>, changed: readonly string[]) =>
  `Changed the ${changed.join(' and ')} of rep ${cmd.index} of ${cmd.object}.`;
