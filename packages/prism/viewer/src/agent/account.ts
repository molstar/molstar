/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SceneNode } from '../model/scene.js';
import type { PrismViewer } from '../viewer.js';
import type { AgentAccountModule } from './tools.js';

// What the agent sees each turn (agent.md §5): the scene's outline (structures with their components and
// representations, volumes, measurements), capped, and the focus (what is selected and focused now: never stored). Text
// that comes from data (a structure's title, an entity's name, a file name) is quoted as data (agent.md §4): in quotes,
// escaped (controls, line and paragraph separators, bidi and zero-width marks gone), capped in length, so it is never
// read as instructions. The person's notes are quoted whole (the viewer's words name data).

/** Lines in the per-turn state block; the rest through `get_scene` (Stories' state block: 40). */
export const StateLines = 24;

/**
 * How much of the scene an outline says: `summary` the structures, volumes and measurements with their components'
 * names; `components` (each turn's) their components too, with their representations and atoms in a line each, and what
 * is focused; `full` everything the tree has (get_scene's default): representations, hidden rows (the unit cell,
 * streaming).
 */
export type OutlineDetail = 'summary' | 'components' | 'full';

/** A label from data, as the model reads it: quoted, escaped, at most 80 characters. */
export function quoted(text: string): string {
  const clean = text
    .replace(/[\u0000-\u001f\u007f-\u009f\u2028\u2029]/g, ' ')
    .replace(/[\u200b-\u200f\u202a-\u202e\u2060-\u2069\ufeff]/g, '')
    .replace(/[`"\\]/g, (c) => `\\${c}`);
  return `"${clean.length > 80 ? `${clean.slice(0, 79)}…` : clean}"`;
}

/** The scene in words, read from the scene's tree (DOM-free) and the interaction sources. */
export function createAccount(viewer: Pick<PrismViewer, 'views' | 'sources'>): AgentAccountModule {
  return {
    // The scene's tree is the Structure panel's model (its feature); without it, no outline.
    outline: (max = StateLines, detail: OutlineDetail = 'full') =>
      'scene' in viewer.views ? outline(viewer.views.scene.roots.value, max, detail) : 'The scene is empty.',
    focus: () => {
      const { selection, focus } = viewer.sources.interaction;
      return [
        `Selected: ${selection.value ? quoted(selection.value) : 'nothing'}.`,
        `Focused: ${focus.value ? quoted(focus.value) : 'nothing'}.`,
      ].join(' ');
    },
  };
}

/** A structure's group of rows (Components, Other): not said itself. */
const isGroup = (n: SceneNode) => n.kind === 'components' || n.kind === 'other';
/** Rows with the groups' rows in their place. */
const ungrouped = (list: readonly SceneNode[]): SceneNode[] =>
  list.flatMap((c) => (isGroup(c) ? ungrouped(c.children.value) : [c]));

/** The tree as indented lines, depth first, without the lazy contents (entities and chains, which no read lists yet). */
function outline(roots: readonly SceneNode[], max: number, detail: OutlineDetail): string {
  if (!roots.length) return 'The scene is empty.';
  const lines: string[] = [];
  let more = 0;
  // What a level leaves out: below `components`, the representations and the rows that say nothing now (hidden
  // features' rows, an empty focus); `summary` keeps the top rows only.
  const said = (n: SceneNode, depth: number) => {
    if (n.kind === 'contents') return false;
    if (detail === 'full') return true;
    if (detail === 'summary') return depth === 0;
    if (n.kind === 'representation' || n.kind === 'volume-representation' || n.kind === 'particles-representation')
      return false;
    if ((n.kind === 'extension' || n.kind === 'unit-cell') && n.isHidden.value) return false;
    return !(n.kind === 'focus' && !n.children.value.length);
  };
  const count = (nodes: readonly SceneNode[], depth: number): number =>
    nodes.reduce(
      (n, c) =>
        isGroup(c)
          ? n + count(c.children.value, depth)
          : said(c, depth)
            ? n + 1 + count(c.children.value, depth + 1)
            : n,
      0,
    );
  const walk = (nodes: readonly SceneNode[], depth: number) => {
    for (const n of nodes) {
      // A structure's groups (Components, Other) are the tree's arrangement: their rows are said as the structure's own.
      if (isGroup(n)) {
        walk(n.children.value, depth);
        continue;
      }
      if (!said(n, depth)) continue;
      if (lines.length >= max) {
        // The items left out, with what is under them.
        more += 1 + count(n.children.value, depth + 1);
        continue;
      }
      // A description is made of labels from data too (a measurement's atoms, a representation's name).
      const description = n.description.value;
      const hidden = n.isHidden.value ? ' (hidden)' : '';
      // A summary names a structure's components on its line.
      const parts =
        detail === 'summary' && n.kind === 'structure'
          ? ungrouped(n.children.value)
              .filter((c) => c.kind === 'component')
              .map((c) => quoted(c.label.value))
          : [];
      lines.push(
        `${'  '.repeat(depth)}- ${n.kind} ${quoted(n.label.value)}${description ? `: ${quoted(description)}` : ''}${hidden}${parts.length ? `; components ${parts.join(', ')}` : ''}`,
      );
      walk(n.children.value, depth + 1);
    }
  };
  walk(roots, 0);
  if (more) lines.push(`… ${more} more (get_scene gives all of it).`);
  else if (detail !== 'full') lines.push('(get_scene gives the representations and the hidden rows too.)');
  return lines.join('\n');
}
