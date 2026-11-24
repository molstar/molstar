/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `search` behavior (_spec/README.md "Behaviors", Search; "Search"): the search model on the viewer
// (`viewer.search`) and its commands, run, next, previous, go, clear and save, so a host's box, the palette, a key, an
// agent and a test drive it the same way. Keys only for stepping (`view.keys`, the canvas's scope: a host's box keeps its
// typing); no chrome: the box is the application's.

import { parseChord } from '@molstar/seqx-viewer/bindings';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import type { SearchStrands } from '@molstar/seqx-core/search';
import { SearchModel } from '../search.js';
import type { SeqViewer } from '../viewer.js';
import type { BehaviorContext, BehaviorHandle, BehaviorProvider } from '../viewer-spec.js';

/** What a search command is given. */
export interface SearchContext {
  /** `search.run`: the pattern, the matcher (`search.matcher` when absent), the strands, what it is limited to. */
  readonly pattern?: string;
  readonly matcher?: string;
  readonly strands?: SearchStrands;
  readonly within?: string;
  /** `search.next`, `search.previous`, `search.go`: zoom to the hit (`false`: pan only). */
  readonly focus?: boolean;
  /** `search.go`: the hit's index. */
  readonly index?: number;
  /** `search.save`: the ref's id and label. */
  readonly id?: string;
  readonly label?: string;
}

export interface SearchBehavior extends BehaviorHandle {
  readonly model: SearchModel;
  readonly commands: readonly CommandDefinition<SearchContext>[];
}

export function searchBehavior(ctx: BehaviorContext): SearchBehavior {
  let viewer: SeqViewer | undefined;
  const model = new SearchModel({
    ...ctx,
    matchers: () => viewer?.context.matchers ?? {},
    notify: (text) => viewer?.notify(text),
  });
  const said = (o: { message: string }): CommandResult => ({ message: o.message });
  const hits = () => (model.hits.value.hits.length ? true : 'No hits');
  const key = (id: string) => {
    const k = parseChord(ctx.config.view.keys[id]?.[0] ?? '');
    return k ? { shortcut: { ...k }, scope: 'canvas' as const } : {};
  };
  const commands: CommandDefinition<SearchContext>[] = [
    {
      id: 'search.run',
      title: 'Search for a pattern',
      category: 'Search',
      when: (c) => (c?.pattern ? true : 'No pattern'),
      run: async (_, c): Promise<CommandResult> =>
        said(
          await model.run({
            pattern: c.pattern!,
            ...(c.matcher !== undefined ? { matcher: c.matcher } : {}),
            ...(c.strands !== undefined ? { strands: c.strands } : {}),
            ...(c.within !== undefined ? { within: c.within } : {}),
          }),
        ),
    },
    {
      id: 'search.next',
      title: 'Next hit',
      category: 'Search',
      ...key('search.next'),
      when: hits,
      run: (_, c): CommandResult => said(model.next({ focus: c?.focus ?? false })),
    },
    {
      id: 'search.previous',
      title: 'Previous hit',
      category: 'Search',
      ...key('search.previous'),
      when: hits,
      run: (_, c): CommandResult => said(model.previous({ focus: c?.focus ?? false })),
    },
    {
      id: 'search.go',
      title: 'Go to a hit',
      category: 'Search',
      when: (c) => (c?.index === undefined ? 'No hit given' : hits()),
      run: (_, c): CommandResult => said(model.go(c.index!, { focus: c.focus ?? false })),
    },
    {
      id: 'search.clear',
      title: 'Clear the search',
      category: 'Search',
      when: () => (model.query.value ? true : 'Nothing searched'),
      run: (): CommandResult => {
        model.clear();
        return { message: 'Search cleared' };
      },
    },
    {
      id: 'search.save',
      title: 'Save the hits as a set',
      category: 'Search',
      when: hits,
      run: async (_, c): Promise<CommandResult> => {
        const o = await model.save({
          ...(c?.id !== undefined ? { id: c.id } : {}),
          ...(c?.label !== undefined ? { label: c.label } : {}),
        });
        return { message: o.message, ...(o.issues.length ? {} : { undoable: true }) };
      },
    },
  ];
  return {
    model,
    commands,
    composed(v) {
      viewer = v;
    },
    dispose() {
      model.retire();
    },
  };
}

/** The search behavior (_spec/README.md "Behaviors", Search): a host lists it, or leaves search out. */
export const SearchProvider: BehaviorProvider<undefined, SearchBehavior> = {
  name: 'search',
  create: (ctx) => searchBehavior(ctx),
};
