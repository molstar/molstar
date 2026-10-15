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

export class PrismViewerSpecError extends Error {
  constructor(readonly problems: readonly string[]) {
    super(`The viewer's spec has problems: ${problems.join('; ')}`);
    this.name = 'PrismViewerSpecError';
  }
}

export function checkSpec(spec: PrismViewerSpec): void {
  const problems: string[] = [];
  const owners = new Set<string>(['core']);
  for (const e of spec.extensions ?? []) {
    if (!IdPattern.test(e.id)) problems.push(`extension id '${e.id}' is not lowercase words joined by dashes`);
    if (owners.has(e.id)) problems.push(`extension '${e.id}' is listed twice`);
    owners.add(e.id);
  }
  for (const id of Object.keys(spec.overrides ?? {})) {
    if (!PartIdPattern.test(id)) {
      problems.push(`override '${id}' is not a part id (<owner>.<kind>.<name>)`);
      continue;
    }
    const owner = id.split('.')[0];
    if (!owners.has(owner)) problems.push(`override '${id}': there is no extension '${owner}'`);
  }
  if (spec.persist && !spec.persist.key) problems.push('persist.key is empty');
  if (problems.length) throw new PrismViewerSpecError(problems);
}
