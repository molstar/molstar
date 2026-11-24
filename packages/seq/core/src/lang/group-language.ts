/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `group` kind (_spec/README.md, Language, Groups of tracks): named groups of a section's tracks, one level deep, folded by
// the reader as session state from the opening `collapsed`. Its own file: `language.ts` is past the module bound.

/** The kinds that hold a section's tracks: the section, and a group of it (one level: a group holds tracks only). */
export const TrackParents = ['section', 'group'] as const;

/**
 * The kinds a reader of an older minor lacks and reads through, each with the version it came in: its children read
 * as its parent's (a `group`'s tracks as its section's, every group open, its `summary` track left out; EC-244).
 */
export const FlattenedKinds: Readonly<Record<string, string>> = { group: '0.11.0' };

export const GroupKinds = {
  group: {
    doc: 'A named group of a section’s tracks (one level: its children are tracks): drawn under a header line, folded by the reader to one summary row.',
    parents: ['section'],
    ids: 'required',
    params: {
      label: { type: 'string', doc: 'The group’s label (its id when absent).' },
      collapsed: {
        type: 'boolean',
        default: false,
        doc: 'Whether the group opens folded: the opening state (session state after: `groups.keep` writes the folds shown here).',
      },
      hidden: {
        type: 'boolean',
        default: false,
        doc: 'Whether the group is hidden: kept in the spec, it and its tracks not drawn; its tracks keep their own `hidden`.',
      },
      summary: {
        type: 'union',
        of: [
          { type: 'string', enum: ['auto', 'none'] },
          { type: 'ref', kinds: ['track'], in: 'group' },
        ],
        doc: 'What the group draws folded: `auto` (when absent: its tracks’ layers in their summary forms, one band), `none` (the header alone), or the id of one of its tracks (drawn only while folded).',
      },
    },
  },
} as const;
