/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Part, PartId } from '../spec.js';

// The spec's overrides (viewer-spec.md §4), one rule for every part, the core's included: `null` removes the part, a
// value replaces its fields (its id kept). Commands go by their command id; views by part id, `<owner>.<kind>.<name>`.
// The core's views keep short ids in the shell (`home`, `reset`), which are `core.<kind>.<id>` as parts.

/** The core's view parts, as overrides and `ui.layout` name them (checked by `spec-check.ts`). */
export const CoreParts: Readonly<Record<string, readonly string[]>> = {
  tab: ['home', 'structure', 'log', 'help', 'sequence'],
  toolbar: ['theme', 'palette', 'reset', 'expand', 'sequence'],
  overlay: ['sequence-legend'],
};

/** The commands the page's shell registers (no headless viewer has them; overrides of them are still valid). */
export const ShellCommandIds: readonly string[] = [
  'help.show',
  'view.theme',
  'sequence.none',
  'sequence.text',
  'sequence.tracks',
  'sequence.rail.ruler',
  'sequence.rail.strip',
  'sequence.rail.band',
  'sequence.rail.none',
  'sequence.legend',
];

/** A core part's id in the shell from its part id (`core.tab.home` → `home`); other ids are the shell's as they are. */
export function shellIdOf(id: string): string {
  if (!id.startsWith('core.')) return id;
  const rest = id.slice('core.'.length);
  const dot = rest.indexOf('.');
  return dot < 0 ? id : rest.slice(dot + 1);
}

export class Overrides {
  readonly #map: Readonly<Record<PartId, Part | null>>;
  readonly #used = new Set<string>();

  constructor(map: Readonly<Record<PartId, Part | null>> | undefined) {
    this.#map = map ?? {};
  }

  get keys(): readonly string[] {
    return Object.keys(this.#map);
  }

  /** Commands by id. */
  commands<T extends { readonly id: string }>(list: readonly T[]): T[] {
    return this.#apply(list, (c) => c.id);
  }

  /** An extension's views, by their part ids. */
  parts<T extends { readonly id: string }>(list: readonly T[] | undefined): T[] {
    return this.#apply(list ?? [], (p) => p.id);
  }

  /** The core's views of a kind, by `core.<kind>.<id>`. */
  core<T extends { readonly id: string }>(kind: string, list: readonly T[]): T[] {
    return this.#apply(list, (p) => `core.${kind}.${p.id}`);
  }

  /** Whether an override removes this part or command. */
  removes(key: string): boolean {
    return this.#map[key] === null;
  }

  /** Override keys nothing has asked for yet. */
  unused(): string[] {
    return this.keys.filter((k) => !this.#used.has(k));
  }

  #apply<T extends { readonly id: string }>(list: readonly T[], keyOf: (x: T) => string): T[] {
    return list.flatMap((x) => {
      const key = keyOf(x);
      if (!(key in this.#map)) return [x];
      this.#used.add(key);
      const replacement = this.#map[key];
      return replacement === null ? [] : [{ ...x, ...(replacement as object), id: x.id } as T];
    });
  }
}
