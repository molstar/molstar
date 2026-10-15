/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PrismViewerSpec } from './spec.js';

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

export function checkSpec(spec: PrismViewerSpec, coreAreas: readonly string[] = []): void {
  const problems: string[] = [];
  const owners = new Set<string>(['core']);
  for (const e of spec.extensions ?? []) {
    if (!IdPattern.test(e.id)) problems.push(`extension id '${e.id}' is not lowercase words joined by dashes`);
    if (owners.has(e.id)) problems.push(`extension '${e.id}' is listed twice`);
    if (coreAreas.includes(e.id)) problems.push(`extension id '${e.id}' is a core API area's`);
    owners.add(e.id);
  }
  for (const id of Object.keys(spec.overrides ?? {})) {
    const parts = id.split('.');
    const isPart = PartIdPattern.test(id) && PartKinds.has(parts[1]);
    if (!isPart && !CommandIdPattern.test(id)) {
      problems.push(`override '${id}' is neither a part id (<owner>.<kind>.<name>) nor a command id`);
      continue;
    }
    if (isPart && !owners.has(parts[0])) problems.push(`override '${id}': there is no extension '${parts[0]}'`);
  }
  if (spec.persist && !spec.persist.key) problems.push('persist.key is empty');
  if (problems.length) throw new PrismViewerSpecError(problems);
}
