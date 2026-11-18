/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `addressing` behavior (_spec/README.md "Behaviors"): the actions on addresses as commands (focus, reveal, select,
// highlight on a named channel, hover on the links' own channel, copy a link, add a ref), so a link, a key, a menu, an
// agent and a test drive the view the same way. A link runs its actions in order, each through the model method its
// command runs.

import type { Issue } from '@molstar/core/repr/tree';
import type { RefInput } from '@molstar/seqx-core/lang';
import type { LinkActionName } from '@molstar/seqx-core/query';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import type { AddRefResult, AddressModel, AddressOutcome } from '@molstar/seqx-viewer/models';
import type { EditContext } from '@molstar/seqx-viewer/ports';

/** What a command on an address is given: the address, and for some a channel, a section or a ref to add. */
export interface AddressContext {
  readonly address?: string;
  /** The highlight channel (`highlight`): the links' channel when absent. */
  readonly channel?: string;
  /** The section a focus moves (each leading section showing it when absent). */
  readonly section?: string;
  /** The ref `refs.add` adds. */
  readonly ref?: RefInput;
  /** The host's run `refs.add` is part of (a host's step that runs it). */
  readonly within?: EditContext;
}

/** The commands' ids by the link action they run. */
export const AddressCommandIds: Readonly<Record<LinkActionName, string>> = {
  focus: 'address.focus',
  reveal: 'address.reveal',
  select: 'address.select',
  highlight: 'address.highlight',
  hover: 'address.hover',
};

export interface LinkRun {
  /** Each action's outcome, in the link's order. */
  readonly outcomes: readonly AddressOutcome[];
  /** The link's own issues (an unknown action, a malformed entry) and each action's. */
  readonly issues: readonly Issue[];
}

export interface AddressingBehavior {
  readonly commands: readonly CommandDefinition<AddressContext>[];
  /** Runs a link's actions through the commands; undefined when the text is no link. */
  runLink(text: string): LinkRun | undefined;
  dispose(): void;
}

export function addressingBehavior(addresses: AddressModel): AddressingBehavior {
  let disposed = false;
  const needsAddress = (ctx: AddressContext): true | string =>
    disposed ? 'The viewer is retired' : ctx.address === undefined ? 'No address' : true;
  const result = (o: AddressOutcome): CommandResult => ({ message: o.message });

  const action = (
    id: string,
    title: string,
    run: (address: string, ctx: AddressContext) => AddressOutcome,
  ): CommandDefinition<AddressContext> => ({
    id,
    title,
    category: 'Address',
    when: needsAddress,
    run: (_, ctx): CommandResult => result(run(ctx.address!, ctx)),
  });

  const commands: CommandDefinition<AddressContext>[] = [
    action(AddressCommandIds.focus, 'Focus an address', (a, ctx) =>
      addresses.focus(a, ctx.section !== undefined ? { section: ctx.section } : {}),
    ),
    action(AddressCommandIds.reveal, 'Reveal an address', (a, ctx) =>
      addresses.reveal(a, ctx.section !== undefined ? { section: ctx.section } : {}),
    ),
    action(AddressCommandIds.select, 'Select an address', (a) => addresses.select(a)),
    action(AddressCommandIds.highlight, 'Highlight an address', (a, ctx) => addresses.highlightAddress(a, ctx.channel)),
    action(AddressCommandIds.hover, 'Hover an address', (a) => addresses.hoverAddress(a)),
    {
      id: 'address.unhover',
      title: 'End a link’s hover',
      category: 'Address',
      run: (): CommandResult => {
        addresses.unhover();
        return { message: '' };
      },
    },
    {
      id: 'address.copy-link',
      title: 'Copy a link',
      category: 'Address',
      // the link is the result: the mounted UI writes it to the clipboard
      when: () => (disposed ? 'The viewer is retired' : addresses.copyLink() !== undefined ? true : 'Nothing to link'),
      run: (): CommandResult => ({ message: addresses.copyLink() ?? '' }),
    },
    {
      id: 'refs.add',
      title: 'Add a bookmark',
      category: 'Address',
      when: (ctx) => (disposed ? 'The viewer is retired' : ctx.ref ? true : 'No ref to add'),
      run: (_, ctx): CommandResult | Promise<CommandResult> => {
        const said = (r: AddRefResult): CommandResult =>
          r === true ? { message: 'Bookmarked', undoable: true } : { message: r[0]?.message ?? 'Refused' };
        const r = addresses.addRef(ctx.ref!, ctx.within ? { within: ctx.within } : {});
        // one step of the host's history when one is bound: it runs in turn
        return r instanceof Promise ? r.then(said) : said(r);
      },
    },
  ];

  return {
    commands,
    runLink(text) {
      if (disposed) return undefined;
      const parsed = addresses.parseLink(text);
      if (!parsed) return undefined;
      const outcomes: AddressOutcome[] = [];
      const issues: Issue[] = [...parsed.issues];
      for (const a of parsed.link?.actions ?? []) {
        const o = addresses.run(a.action, a.address, a.action === 'highlight' ? parsed.link?.channel : undefined);
        outcomes.push(o);
        issues.push(...o.issues);
      }
      return { outcomes, issues };
    },
    dispose() {
      disposed = true;
    },
  };
}
