/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { AntibodyViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { complexEditStory } from '../complex/edit-root.js';
import { exampleStory } from '../example-root.js';
import { bookmark, bookmarkLines, readoutLines, recolor, remove, rename, selectCdr, setRuler } from './compose.js';

const view = (id: string) => AntibodyViews.find((v) => v.id === id)!;

const meta: Meta = { title: 'Examples/Antibody/Views' };
export default meta;

/**
 * Trastuzumab's light and heavy chains (PDB 1N8Z), a Kabat ruler and an IMGT ruler over each, the CDRs by both
 * definitions below: the heavy chain's Kabat insertions 52A, 82A–C and 100A–C. Hover a residue: the viewer's card says
 * its label in both numberings (the panel, once the hover's readout, names the selection alone). Editing: switch which scheme the top ruler shows (a patch of its layer's `numbering`).
 */
export const Schemes: StoryObj = {
  render: () =>
    complexEditStory({
      view: view('schemes'),
      lines: readoutLines,
      edits: [
        { label: 'Top rulers in IMGT', run: (s) => setRuler(s, 'imgt') },
        { label: 'Top rulers in Kabat', run: (s) => setRuler(s, 'kabat') },
        { label: 'Top rulers in ordinals', run: (s) => setRuler(s, 'ordinal') },
      ],
    }),
};

/** The design under the same rulers: IMGT's 111.1 and 112.1 in CDR-H3, before 112 (EC-7); Kabat's 100D and 100E. */
export const DesignNumberings: StoryObj = {
  render: () => complexEditStory({ view: view('schemes-design'), lines: readoutLines, edits: [] }),
};

/** Each chain's design and wild type as two rows by their IMGT positions: the insertions and the deletion as gaps. */
export const DesignAgainstWildType: StoryObj = { render: () => exampleStory(view('against')) };

/**
 * The CDR bookmarks: select CDR-H3 (or L1) by its Kabat address, bookmark the selection, rename it, recolor it, remove
 * it; each one step, Undo and Redo walk them. The panel lists the bookmarks.
 */
export const CdrBookmarks: StoryObj = {
  render: () =>
    complexEditStory({
      view: view('cdrs'),
      lines: bookmarkLines,
      edits: [
        { label: 'Select CDR-H3', run: (s) => selectCdr(s, 'heavy', 3) },
        { label: 'Select CDR-L1', run: (s) => selectCdr(s, 'light', 1) },
        { label: 'Bookmark the selection', run: (s) => bookmark(s, 'cdr', 'CDR') },
        { label: 'Rename it', run: (s) => rename(s, 'cdr', 'The CDR, renamed') },
        { label: 'Recolor it', run: (s) => recolor(s, 'cdr', '#d62728') },
        { label: 'Remove it', run: (s) => remove(s, 'cdr') },
      ],
    }),
};
