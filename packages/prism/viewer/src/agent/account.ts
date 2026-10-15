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
// that comes from data (a structure's title, an entity's name, a file name) is quoted as data (agent.md §4): in a
// fence, escaped, capped in length, so it is never read as instructions.

/** Lines in the per-turn state block; the rest through `get_scene` (Stories' state block: 40). */
export const StateLines = 40;

/** A label from data, as the model reads it: quoted, escaped, at most 80 characters. */
export function quoted(text: string): string {
  const clean = text.replace(/[\u0000-\u001f\u007f]/g, ' ').replace(/[`"\\]/g, (c) => `\\${c}`);
  return `"${clean.length > 80 ? `${clean.slice(0, 79)}…` : clean}"`;
}

/** The scene in words, read from the scene's tree (DOM-free) and the interaction sources. */
export function createAccount(viewer: Pick<PrismViewer, 'views' | 'sources'>): AgentAccountModule {
  return {
    outline: (max = StateLines) => outline(viewer.views.scene.roots.value, max),
    focus: () => {
      const { selection, focus } = viewer.sources.interaction;
      return [
        `Selected: ${selection.value ? quoted(selection.value) : 'nothing'}.`,
        `Focused: ${focus.value ? quoted(focus.value) : 'nothing'}.`,
      ].join(' ');
    },
  };
}

/** The tree as indented lines, depth first, without the lazy contents (entities and chains: `get_scene` reads them). */
function outline(roots: readonly SceneNode[], max: number): string {
  if (!roots.length) return 'The scene is empty.';
  const lines: string[] = [];
  let more = 0;
  const walk = (nodes: readonly SceneNode[], depth: number) => {
    for (const n of nodes) {
      if (n.kind === 'contents') continue;
      if (lines.length >= max) {
        more++;
        continue;
      }
      const description = n.description.value;
      const hidden = n.isHidden.value ? ' (hidden)' : '';
      lines.push(
        `${'  '.repeat(depth)}- ${n.kind} ${quoted(n.label.value)}${description ? `: ${description}` : ''}${hidden}`,
      );
      walk(n.children.value, depth + 1);
    }
  };
  walk(roots, 0);
  if (more) lines.push(`… ${more} more (get_scene gives all of it).`);
  return lines.join('\n');
}
