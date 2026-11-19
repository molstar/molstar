/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type ReadonlyState, State } from '@molstar/core/reactive';
import { formatChord } from '@molstar/ui/commands/keys';
import { parseChord } from './bindings.js';
import type { SeqViewer } from './viewer.js';

// The default toolbar's buttons (_spec/README.md "Input", Toolbar): DOM-free, so the headless twin reads what the
// mounted bar shows. A button runs a command of the viewer by id; its title says what it does with its key
// (`view.keys`, the first chord); it is disabled, with the reason, while its command cannot run. Whether a command can
// run is read from the models the command reads and is reactive for the view's history (`view.history`): the one pair
// of buttons the bar has by default; a command with no reactive reason is never disabled.

/** One button of the bar. */
export interface ToolbarButton {
  /** The command it runs. */
  readonly command: string;
  /** An icon of `@molstar/ui`'s set. */
  readonly icon: string;
  /** What it does, and its key: the accessible name is `label`, the tooltip is this. */
  readonly title: string;
  readonly label: string;
  /** The first chord of the command's keys as shown on this platform (`Alt+←`), when it has one. */
  readonly key: string | undefined;
  /** Why it cannot run now, else false. */
  readonly disabled: ReadonlyState<string | false>;
  /** Runs the command on the viewer's first leading section; its message. */
  run(): string | Promise<string>;
}

/** The commands the bar can show, with their icons and labels. */
const Known: Readonly<Record<string, { readonly icon: string; readonly label: string }>> = {
  'view.back': { icon: 'chevron-left', label: 'Back' },
  'view.forward': { icon: 'chevron-right', label: 'Forward' },
  // where to find what the view answers to, said once, quietly: its title is "Keys and gestures (?)"
  'help.keys': { icon: 'keyboard', label: 'Keys and gestures' },
};

/** The commands the default bar shows (`SeqViewerSpec.ui.toolbar`). */
export const DefaultToolbar: readonly string[] = ['view.back', 'view.forward', 'help.keys'];

/** Never disabled: a command with no reason not to run (the keys panel's). */
const Enabled: ReadonlyState<string | false> = new State<string | false>(false).asReadonly;

/** The buttons of a bar for the commands named (unknown ones are left out), over a viewer's models. */
export function toolbarButtons(viewer: SeqViewer, commands: readonly string[] = DefaultToolbar): ToolbarButton[] {
  const out: ToolbarButton[] = [];
  for (const command of commands) {
    const known = Known[command];
    // a command the viewer lacks (the keys panel not listed) has no button
    if (!known || !viewer.commands.some((c) => c.id === command)) continue;
    const chord = parseChord(viewer.config.view.keys[command]?.[0] ?? '');
    const key = chord ? formatChord(chord) : undefined;
    const history = viewer.view.history;
    const disabled: ReadonlyState<string | false> =
      command === 'view.back'
        ? derive(history, (h) => (h.back > 0 ? false : 'No view to go back to'))
        : command === 'view.forward'
          ? derive(history, (h) => (h.forward > 0 ? false : 'No view to go forward to'))
          : Enabled;
    out.push({
      command,
      icon: known.icon,
      label: known.label,
      key,
      title: key ? `${known.label} (${key})` : known.label,
      disabled,
      run: () => viewer.execute(command, { section: viewer.view.resolved.value?.sections[0]?.address }),
    });
  }
  return out;
}
