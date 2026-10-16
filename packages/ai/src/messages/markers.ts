/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The words and labels the session writes into a conversation (tools-and-session.md §4), here so the trace reads them
// without importing the session (and, through it, the tools).

/** Opens the person's words in a request: compaction and the trace find them by it. */
export const REQUEST_MARKER = 'Request: ';
/** Opens the Focus line of a request. */
export const FOCUS_MARKER = 'Focus: ';
/** Heads guidance sent again because it changed since the conversation began. */
export const GUIDANCE_UPDATE = 'Updated guidance (it replaces the earlier version):';
/** Heads the host's notes in a request (the default; a host sets its own). */
export const NOTES_HEADING = 'Notes since your last turn:';

/**
 * The labels the session gives the entries it writes (an `EntryLabel`'s `kind`): what compaction, `catchUp` and a
 * trace find them by. The markers above stay in the text, for the model.
 */
/** A turn's request; its data is `{ prompt, focus? }`. */
export const LABEL_REQUEST = 'session:request';
/** The guidance; its key is the text's fingerprint. */
export const LABEL_GUIDANCE = 'session:guidance';
/** Guidance sent again because it changed (`GUIDANCE_UPDATE` and the text); its key is the text's fingerprint. */
export const LABEL_GUIDANCE_UPDATE = 'session:guidance-update';
/** What compaction carried over from the earlier conversation. */
export const LABEL_CARRY_OVER = 'session:carry-over';
/** The host's memory, as compaction placed it. */
export const LABEL_MEMORY = 'session:memory';
