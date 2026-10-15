/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CoreParts } from './extensions/overrides.js';
import { isLazyExtension, type PrismViewerSpec } from './spec.js';

// The first step of creation (viewer-spec.md §6): the spec checked before anything is made. Every problem at once, in
// one error naming each.

const IdPattern = /^[a-z][a-z0-9-]*$/;
const PartIdPattern = /^[a-z][a-z0-9-]*\.[a-z][a-z0-9-]*\.[a-z0-9][a-z0-9.-]*$/;
/** A command's id (`<area>.<verb>`, viewer-spec.md §4: command ids are their own). */
const CommandIdPattern = /^[a-z][a-z0-9-]*(\.[a-z0-9][a-z0-9-]*)+$/;
const PartKinds = new Set([
  'tab',
  'section',
  'row',
  'toolbar',
  'overlay',
  'command',
  'tool',
  'skill',
  'sequence',
  'api',
  'settings',
]);

export class PrismViewerSpecError extends Error {
  constructor(readonly problems: readonly string[]) {
    super(`The viewer's spec has problems: ${problems.join('; ')}`);
    this.name = 'PrismViewerSpecError';
  }
}

/** A part id (`<owner>.<kind>.<name>`), as opposed to a command id. */
export function isPartId(id: string): boolean {
  return PartIdPattern.test(id) && PartKinds.has(id.split('.')[1]);
}

/** Area ids the core has besides `CoreAreas`: made when the spec asks for them. */
const OptionalCoreAreas = ['sequence'];

export function checkSpec(spec: PrismViewerSpec, coreAreas: readonly string[] = []): void {
  const problems: string[] = [];
  const owners = new Set<string>(['core']);
  // What a lazy extension declares; an eager one's views are known when they register (the host checks those).
  const declared = new Map<string, Set<string> | undefined>();
  for (const e of spec.extensions ?? []) {
    if (!IdPattern.test(e.id)) problems.push(`extension id '${e.id}' is not lowercase words joined by dashes`);
    if (owners.has(e.id)) problems.push(`extension '${e.id}' is listed twice`);
    if (coreAreas.includes(e.id) || OptionalCoreAreas.includes(e.id)) {
      problems.push(`extension id '${e.id}' is a core API area's`);
    }
    owners.add(e.id);
    declared.set(e.id, isLazyExtension(e) ? new Set(e.stub.parts.map((p) => p.id)) : undefined);
  }
  for (const id of Object.keys(spec.overrides ?? {})) {
    const isPart = isPartId(id);
    if (!isPart && !CommandIdPattern.test(id)) {
      problems.push(`override '${id}' is neither a part id (<owner>.<kind>.<name>) nor a command id`);
      continue;
    }
    if (!isPart) continue;
    const [owner, kind] = id.split('.');
    const name = id.slice(owner.length + kind.length + 2);
    if (!owners.has(owner)) problems.push(`override '${id}': there is no extension '${owner}'`);
    else if (owner === 'core' && !CoreParts[kind]?.includes(name))
      problems.push(`override '${id}': the core has no such part`);
    else if (declared.get(owner) && !declared.get(owner)!.has(id)) {
      problems.push(`override '${id}': '${owner}' declares no such part`);
    }
  }
  if (spec.persist && !spec.persist.key) problems.push('persist.key is empty');
  if (problems.length) throw new PrismViewerSpecError(problems);
}
