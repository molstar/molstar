/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { formatChord } from '@molstar/prism-ui/commands/keys';

// Help's prose (parity.md §5, "Help (prose, bindings)"): the basics in a few lines, above the bindings in the Help tab
// (`HelpPanel`'s `intro`). Short on purpose: the command menu is where everything else is found.

const Basics: readonly { title: string; text: () => string }[] = [
  {
    title: 'Load',
    text: () =>
      'Drop files on the viewer, or use Home: a PDB or AlphaFold DB entry, a URL, Open files. A session file replaces what is loaded.',
  },
  {
    title: 'Select',
    text: () =>
      'Switch the canvas to Select and click; ⇧ adds. Or run Select… with a target: A:120-130, ligand, resn HEM. The selection bar acts on it.',
  },
  {
    title: 'Focus',
    text: () =>
      'In Interact, click a residue to focus it with its surroundings; click empty space to leave. Double-click empty space resets the camera.',
  },
  {
    title: 'Commands',
    text: () =>
      `${formatChord({ key: 'k', mod: true })} opens the command menu: every action by name, with its keys. ? opens this Help.`,
  },
  {
    title: 'Undo',
    text: () =>
      `${formatChord({ key: 'z', mod: true })} undoes a step, ${formatChord({ key: 'z', mod: true, shift: true })} redoes it; History lists the steps.`,
  },
  {
    title: 'Assistant',
    text: () =>
      'Where the viewer has one, ask in plain words ("color chain A by pLDDT"). It acts through the same commands, each step one Undo.',
  },
];

const Caption =
  'margin: 0; padding: 0 var(--msp-inset) var(--msp-space-1); color: var(--msp-color-text-subtle); font-size: var(--msp-font-size-xs); font-weight: 500';
const Row = 'display: grid; padding: var(--msp-space-1) var(--msp-inset)';

/** The basics: load, select, focus, the command menu, undo, the assistant. */
export function HelpProse() {
  return (
    <section style="display: grid; gap: var(--msp-space-1)" aria-label="Basics">
      <h3 style={Caption}>Basics</h3>
      <dl style="display: grid; margin: 0">
        {Basics.map((b) => (
          <div style={Row}>
            <dt style="color: var(--msp-color-text); font-weight: 600">{b.title}</dt>
            <dd style="margin: 0; color: var(--msp-color-text-muted); line-height: 1.4">{b.text()}</dd>
          </div>
        ))}
      </dl>
    </section>
  );
}
