/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Issue } from '@molstar/core/repr/tree';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { SeqViewer } from './viewer.js';

// The model's issues as the viewer says them (_spec/README.md "SpecModel", Issues of a pending address): an issue of a
// pending address is not shown. `seq.resolve.length.unknown` of a system read from a source is transient while the
// source opens (its coverage gives the length, the model is resolved again and the issue goes), so it is held back
// until the source's opening settles: opened (the length still unknown: then it is said) or failed (its own issue).
// Every other issue is shown as the model has it; `seq.resolve.system.unnamed` (the source opened and names no system
// after the sequence) is never transient.

/** The code held back while its system's source is pending. */
const Transient = 'seq.resolve.length.unknown';

/**
 * The model's issues less those of a pending address: `seq.resolve.length.unknown` at an annotation or a map over a
 * system whose source `pending` says has not settled.
 */
export function shownModelIssues(model: StaticModel, pending: (source: string) => boolean): readonly Issue[] {
  if (!model.issues.some((i) => i.code === Transient)) return model.issues;
  return model.issues.filter((i) => i.code !== Transient || !waits(model, i.address, pending));
}

/** Whether a system the issue's node reads (an annotation's target, a map's ends) is read from a pending source. */
function waits(model: StaticModel, address: string | undefined, pending: (source: string) => boolean): boolean {
  if (address === undefined) return false;
  const map = model.maps.get(address);
  const systems = map ? [map.from, map.to] : [model.annotations.get(address)?.target];
  for (const s of systems) {
    const source = s === undefined ? undefined : model.systems.get(s)?.source;
    if (source !== undefined && pending(source)) return true;
  }
  return false;
}

/**
 * The viewer's model's issues as it says them: a source is pending until the viewer's store opened it from the model's
 * own entry, or reported it (a failed opening).
 */
export function modelIssues(viewer: SeqViewer): readonly Issue[] {
  const model = viewer.spec.model.value;
  const reported = viewer.store.issues.value;
  return shownModelIssues(
    model,
    (source) => viewer.store.openedCoverage(model, source) === undefined && !reported.some((i) => i.address === source),
  );
}
