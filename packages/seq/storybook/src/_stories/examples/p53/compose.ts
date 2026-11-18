/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Patch } from '@molstar/core/repr/patch';
import type { Node } from '@molstar/core/repr/tree';
import { printLink } from '@molstar/seqx-core/query';
import type { LinkRun } from '@molstar/seqx-viewer/behaviors';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { followLink } from '../../features/addressing/compose.js';
import { type EditingSession, proteinSection, SectionPath, track } from '../../features/editing/compose.js';

// The p53 example stories' shared composition, DOM-free: each story and its headless twin call these. The pieces are the
// Features stories' (editing, addressing, selection); what is here is what the examples add on p53. Every action goes
// through the viewer's public API (commands, `spec.patch`, its input), never around it.

/** The id of the track the example adds. */
export const AddedTrack = 'domains-2';

/**
 * Adds a second domains track after the first, a copy of it with its own id and label: one patch, so one step of the
 * host's stack. A second call finds the id taken and says so.
 */
export async function addTrack(s: EditingSession): Promise<string> {
  const section = proteinSection(s);
  const domains = section?.children?.find((c) => c.id === 'domains');
  if (!section || !domains) return 'No domains track to copy';
  if (section.children?.some((c) => c.id === AddedTrack)) return 'Already added';
  const node: Node = { ...domains, id: AddedTrack, params: { ...domains.params, label: 'Domains (added)' } };
  const patch: Patch = { ops: [{ op: 'insert', parent: SectionPath, place: { after: 'domains' }, node }] };
  const r = await s.viewer.spec.patch('Add a track', patch);
  return r.ok ? 'Added' : `Not added: ${r.issues.map((i) => i.message).join('; ')}`;
}

/** The track ids of the protein section, in order. */
export function trackOrder(s: EditingSession): string[] {
  return (proteinSection(s)?.children ?? []).map((c) => c.id as string);
}

/** A track's height as the spec shown says (a preview included), else as laid out; undefined before a layout. */
export function trackHeight(s: EditingSession, id: string): number | undefined {
  const set = proteinSection(s, true)?.children?.find((c) => c.id === id)?.params?.height;
  if (typeof set === 'number') return set;
  const laid = s.viewer.frame?.layout.sections.flatMap((x) => x.tracks).find((t) => t.track.address === track(id));
  return laid?.box.height;
}

/** What a link run says: its outcomes' messages, and the issues' only when no outcome waits (a wait says it once). */
export function describeLinkRun(run: LinkRun | undefined): string {
  if (!run) return 'Not a link';
  const said = run.outcomes.map((o) => o.message).filter(Boolean);
  const waiting = run.outcomes.some((o) => o.status === 'waiting');
  const more = waiting ? [] : run.issues.map((i) => i.message);
  return [...new Set([...said, ...more])].join('; ');
}

/**
 * A link's status line that follows a wait: a link on an address that arrives later says it waits, and `Focused` once
 * the node has arrived and the view has moved to it. A later link, or `dispose`, ends the watch. A link followed before
 * `loaded` (the view loaded and settled) is held and runs then: until a spec is loaded no address names anything, and the viewer answers `Nothing to
 * focus` rather than wait.
 */
export class LinkStatus {
  private _off: (() => void) | undefined;
  private _loaded = false;
  private _held: string | undefined;

  constructor(
    private readonly viewer: SeqViewer,
    private readonly say: (message: string) => void,
  ) {}

  /** The view is loaded: a held link runs now. */
  loaded(): void {
    this._loaded = true;
    const href = this._held;
    this._held = undefined;
    if (href !== undefined) this.say(this.follow(href));
  }

  /** Runs the link: the message now, and the later one through `say`. */
  follow(href: string): string {
    this._stop();
    if (!this._loaded) {
      this._held = href;
      return 'Waiting for the view to load';
    }
    const run = followLink(this.viewer, href);
    if (run?.outcomes.some((o) => o.status === 'waiting')) {
      const address = this.viewer.addresses.parseLink(href)?.link?.actions[0]?.address;
      if (address !== undefined) {
        this._off = this.viewer.rested.listen(() => {
          if (this.viewer.addresses.resolve(address).status !== 'resolved') return;
          if (this.viewer.addresses.focusAddress() !== address) return;
          this._stop();
          this.say('Focused');
        });
      }
    }
    return describeLinkRun(run);
  }

  dispose(): void {
    this._stop();
  }

  private _stop(): void {
    this._off?.();
    this._off = undefined;
  }
}

/** A run of the short text: plain, or a link (its text and its `!seq-…` href). */
export type Run = string | { readonly text: string; readonly href: string };

const link = (text: string, address: string): Run => ({
  text,
  href: printLink({ actions: [{ action: 'focus', address }] }),
});

/** The text beside the gene's region: each link moves the view; a transcript's arrives later, and its focus waits. */
export const RegionText: readonly Run[] = [
  'TP53 sits at ',
  link('the gene', '/refs/gene'),
  ' of this region. Its second transcript, ',
  link('T2', '/refs/t2'),
  ', is read with its part of the region, so a focus on it waits until it arrives; so does one on a neighbor, ',
  link('WRAP53', '/refs/wrap53'),
  '.',
];
