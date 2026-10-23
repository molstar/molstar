/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet } from '@molstar/core/data/int';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type RuntimeContext, Task } from '@molstar/core/task';
import { type Model, Structure, type Trajectory } from '@molstar/model/model/structure';
import type { AtomicHierarchy } from '@molstar/model/model/structure/model/properties/atomic';
import { ModelPeriodicity } from '@molstar/model/model/structure/model/periodic';
import { Script } from '@molstar/model/script/script';

// What a frame in a simulation box shows (model.md §4 "Periodic frames"): `ModelFromTrajectory`'s params.

export type Periodic = 'whole' | 'box' | 'nojump' | 'off';

/** Polymers: the protein (and nucleic acids) of a simulation. */
export const PolymerCenter = Script('(sel.atom.atom-groups :entity-test (= atom.entity-type polymer))', 'mol-script');

export const PeriodicParams = {
  periodic: PD.Optional(
    PD.Select<Periodic>(
      'whole',
      [
        ['whole', 'Molecules whole'],
        ['box', 'Molecules in the box'],
        ['nojump', 'No jumps between frames'],
        ['off', 'As in the file'],
      ],
      {
        label: 'Periodic box',
        description:
          'For a frame in a simulation box: each molecule split across the box made whole again; also moved into the box; or kept from jumping a box length between frames (molecules leave the box as they diffuse); or the coordinates as in the file',
      },
    ),
  ),
  periodicCenter: PD.Optional(
    PD.MappedStatic(
      'box',
      {
        box: PD.EmptyGroup({ label: 'The box' }),
        selection: PD.Group(
          { script: PD.Script(PolymerCenter) },
          { label: 'A selection', description: 'At the box center, every other molecule around it' },
        ),
      },
      {
        label: 'Center on',
        description: "What the box is centered on: the box as it is, or a selection's atoms (polymers by default)",
        hideIf: (p: { periodic?: Periodic }) => p.periodic !== 'box',
      },
    ),
  ),
};

export type PeriodicProps = {
  periodic?: Periodic;
  periodicCenter?: { name: 'box'; params: Record<string, never> } | { name: 'selection'; params: { script: Script } };
};

// a selection's atoms, once per topology (its hierarchy) and script: their identity keys the frames made
const centers = new WeakMap<AtomicHierarchy, Map<string, ModelPeriodicity.Center | null>>();

/** The selection's atoms in the frame's topology, evaluated once on its first frame (made whole); none: undefined. */
function centerOf(frame: Model, script: Script): ModelPeriodicity.Center | undefined {
  let byScript = centers.get(frame.atomicHierarchy);
  if (!byScript) {
    byScript = new Map();
    centers.set(frame.atomicHierarchy, byScript);
  }
  const key = `${script.language}:${script.expression}`;
  if (!byScript.has(key)) {
    const loci = Script.toLoci(script, Structure.ofModel(ModelPeriodicity.makeWhole(frame)));
    const atoms: number[] = [];
    for (const e of loci.elements) OrderedSet.forEach(e.indices, (j) => atoms.push(e.unit.elements[j]));
    atoms.sort((a, b) => a - b);
    byScript.set(key, atoms.length ? { atoms } : null);
  }
  return byScript.get(key) ?? undefined;
}

/**
 * The trajectory's frame as the params show it: absent (an older state) is whole, the default. No jump walks the
 * frames before it (`ModelPeriodicity.applyNoJump`).
 */
export async function periodicFrame(
  trajectory: Trajectory,
  index: number,
  props: PeriodicProps,
  ctx: RuntimeContext,
): Promise<Model> {
  if (props.periodic === 'nojump') return ModelPeriodicity.applyNoJump(trajectory, index, ctx);
  return applyPeriodic(await Task.resolveInContext(trajectory.getFrameAtIndex(index), ctx), props);
}

/** The frame as the params show it (but no jump, which needs the trajectory: `periodicFrame`). */
export function applyPeriodic(frame: Model, props: PeriodicProps): Model {
  const periodic = props.periodic ?? 'whole';
  if (periodic === 'off' || periodic === 'nojump' || !ModelPeriodicity.isApplicable(frame)) return frame;
  if (periodic === 'whole') return ModelPeriodicity.makeWhole(frame);
  const c = props.periodicCenter;
  const center = c?.name === 'selection' ? centerOf(frame, c.params.script) : undefined;
  return ModelPeriodicity.apply(frame, center ? { kind: 'box', center } : ModelPeriodicity.InBox);
}

/**
 * The trajectory with every frame as the params show it (`StructureFromTrajectory`: an ensemble of all its frames); the
 * trajectory itself when off or without a simulation box.
 */
export function periodicTrajectory(trajectory: Trajectory, props: PeriodicProps): Trajectory {
  if (props.periodic === 'off' || !ModelPeriodicity.isApplicable(trajectory.representative)) return trajectory;
  const representative =
    props.periodic === 'nojump'
      ? ModelPeriodicity.makeWhole(trajectory.representative)
      : applyPeriodic(trajectory.representative, props);
  return {
    duration: trajectory.duration,
    frameCount: trajectory.frameCount,
    representative,
    getFrameAtIndex: (i) => Task.create('Periodic frame', (ctx) => periodicFrame(trajectory, i, props, ctx)),
  };
}

/** The periodic params, hidden when the trajectory has no simulation box. */
export function periodicParamsFor(trajectory: Trajectory | undefined): typeof PeriodicParams {
  if (!trajectory || ModelPeriodicity.isApplicable(trajectory.representative)) return PeriodicParams;
  return {
    periodic: { ...PeriodicParams.periodic, isHidden: true },
    periodicCenter: { ...PeriodicParams.periodicCenter, isHidden: true },
  };
}
