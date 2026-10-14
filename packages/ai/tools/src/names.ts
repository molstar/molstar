/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** The disclosure tools' default ids (tools-and-session.md §3); a toolkit's `names` override them. */
export const DISCLOSE_SKILLS = 'disclose_skills';
export const DISCLOSE_TOOLS = 'disclose_tools';

/**
 * The labels the toolkit gives the entries it writes (an `EntryLabel`'s `kind`): what a session, compaction and a
 * trace find them by, never by their text.
 */
export const LABEL_DEFINITIONS = 'tools:definitions';
/** A skill's guidance; its key is the group's key (`skill:<id>@<fingerprint>`). */
export const LABEL_SKILL = 'tools:skill';
/** A disclosable entry's text; its key is the entry's key (`entry:<id>@<fingerprint>`). */
export const LABEL_ENTRY = 'tools:entry';
