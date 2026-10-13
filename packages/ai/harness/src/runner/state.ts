/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hashMessage, type Json, type Message } from '@molstar/ai-agent';
import type { ItemView } from '../contracts/drive.js';
import type { RetentionSnapshot } from '../contracts/inspect.js';
import type { Move, RunnerPort } from '../contracts/ports.js';
import type {
  AttemptRecord,
  BlobRef,
  CallId,
  Compaction,
  EntryId,
  HeadId,
  HeadRecord,
  Limits,
  ScopeId,
  Segment,
} from '../contracts/records.js';
import { buildMove } from '../head/index.js';

/** How a step stops early (harness.md §6); the first of stop, fail, park and interrupt wins. */
export class Exit extends Error {
  readonly kind: 'stop' | 'fail' | 'park' | 'interrupt';
  readonly reason?: Json;
  readonly path: readonly CallId[];
  readonly note?: Json;
  readonly item?: ItemView;
  constructor(
    kind: Exit['kind'],
    options: { reason?: Json; path?: readonly CallId[]; note?: Json; item?: ItemView } = {},
  ) {
    super(`exit: ${kind}`);
    this.kind = kind;
    this.path = options.path ?? [];
    if (options.reason !== undefined) this.reason = options.reason;
    if (options.note !== undefined) this.note = options.note;
    if (options.item !== undefined) this.item = options.item;
  }
}
const RANK: Record<Exit['kind'], number> = { stop: 0, fail: 1, park: 2, interrupt: 3 };
export const firstExit = (exits: readonly Exit[]): Exit | undefined =>
  [...exits].sort((a, b) => RANK[a.kind] - RANK[b.kind])[0];

/** Work ends without a commit: the turn is gone (reset), the lease is lost, or the runner was aborted. */
export class Halt extends Error {
  readonly kind: 'reset' | 'lost' | 'abort';
  constructor(kind: Halt['kind']) {
    super(`halt: ${kind}`);
    this.kind = kind;
  }
}

/** What one top-level step's runs share: its abort, its exits, and what an interrupt asked. */
export class StepState {
  readonly controller = new AbortController();
  readonly exits: Exit[] = [];
  /** Open decisions the step must see answered before it settles (an interrupt's). */
  readonly pending: Promise<unknown>[] = [];
  /** The text unstarted calls get instead of running (an interrupt answered `skip`). */
  skip?: string;
  /** The turn's limits when the step began; a commit stores them again if an `extend` replaced them. */
  limits?: Limits;
  /** The step is committing: an interrupt seen now belongs to the next step. */
  closed = false;
  /** Records an exit and aborts the step's other work (harness.md §6). */
  exit(exit: Exit): Exit {
    this.exits.push(exit);
    if (!this.controller.signal.aborted) this.controller.abort(exit);
    return exit;
  }
}

/** A call record while its step is open. */
export interface OpenCall {
  readonly id: CallId;
  readonly key: string;
  readonly parentCallId?: CallId;
  readonly kind: 'model' | 'tool';
  readonly request: BlobRef;
  readonly attempts: AttemptRecord[];
  result?: BlobRef;
}

/** The runner's working set (harness.md §6): immutable records, cached without invalidation. */
export class WorkingSet {
  readonly messages = new Map<EntryId, Message>();
  readonly segments = new Map<string, Segment>();
  readonly heads = new Map<HeadId, HeadRecord>();
  readonly blobs = new Map<BlobRef, Json>();
  private readonly segmentRoots = new Set<string>();
  private readonly entryRefs = new Map<EntryId, number>();
  private readonly tokens = new Map<EntryId, number>();
  /** Messages known to be stored (read, or committed), cached or spilled: a commit does not resend them. */
  readonly stored = new Set<EntryId>();
  /** Evicted records, least recently used first, within `cacheBytes` (keys: `m:`, `s:`, `h:` and the id). */
  private readonly spill = new Map<string, { value: Message | Segment | HeadRecord; size: number }>();
  private spillBytes = 0;
  private readonly port: RunnerPort;
  private readonly cacheBytes: number;
  /** Reads are retried: a failed read is not a loss (harness.md §4). */
  private readonly read: <T>(read: () => Promise<T>) => Promise<T>;

  constructor(port: RunnerPort, cacheBytes: number, read: <T>(read: () => Promise<T>) => Promise<T>) {
    this.port = port;
    this.cacheBytes = cacheBytes;
    this.read = read;
  }

  readSegments = (ids: readonly string[]): Segment[] =>
    ids.map((id) => {
      const segment = this.segments.get(id);
      if (segment === undefined) throw new Error(`segment ${id} is not in the working set`);
      return segment;
    });

  /** Loads a head (record, segments and messages), reading only what is neither cached nor spilled. */
  async load(id: HeadId): Promise<{ record: HeadRecord; entries: EntryId[] }> {
    let record = this.heads.get(id) ?? this.unspill<HeadRecord>(`h:${id}`);
    if (record === undefined) [record] = await this.read(() => this.port.headRecords([id]));
    if (record === undefined) throw new Error(`no head ${id}`);
    this.heads.set(id, record);
    const missing = record.segments.filter((s) => {
      const spilled = this.segments.has(s) ? undefined : this.unspill<Segment>(`s:${s}`);
      if (spilled !== undefined) this.segments.set(s, spilled);
      return !this.segments.has(s);
    });
    if (missing.length > 0)
      for (const segment of await this.read(() => this.port.segments(missing))) this.segments.set(segment.id, segment);
    const entries = this.readSegments(record.segments).flatMap((segment) => segment.entries);
    await this.ensure(entries);
    return { record, entries };
  }

  /** Makes sure the messages are cached: from the spill, else one batched read. */
  async ensure(entries: readonly EntryId[]): Promise<void> {
    const missing = [...new Set(entries)].filter((id) => {
      if (this.messages.has(id)) return false;
      const spilled = this.unspill<Message>(`m:${id}`);
      if (spilled !== undefined) this.messages.set(id, spilled);
      return spilled === undefined;
    });
    if (missing.length === 0) return;
    const found = await this.read(() => this.port.entries(missing));
    missing.forEach((id, i) => {
      const message = found[i];
      if (message !== undefined) this.messages.set(id, message);
      this.stored.add(id);
    });
  }

  /**
   * An input-token estimate of a head (a quarter of its JSON length), cached per entry. A deferred
   * tool declaration counts its name only: a provider shows at most that.
   */
  tokensOf(entries: readonly EntryId[]): number {
    let total = 0;
    for (const id of entries) {
      let tokens = this.tokens.get(id);
      if (tokens === undefined) {
        const message = this.messages.get(id);
        const shown =
          message?.role === 'system' && message.parts.some((p) => p.kind === 'tool-definition' && p.deferred)
            ? message.parts.map((p) => (p.kind === 'tool-definition' && p.deferred ? p.definition.name : p))
            : message;
        tokens = Math.ceil(JSON.stringify(shown ?? null).length / 4);
        this.tokens.set(id, tokens);
      }
      total += tokens;
    }
    return total;
  }

  messagesOf(entries: readonly EntryId[]): Message[] {
    return entries.map((id) => {
      const message = this.messages.get(id);
      if (message === undefined) throw new Error(`message ${id} is not in the working set`);
      return message;
    });
  }

  async add(messages: readonly Message[]): Promise<EntryId[]> {
    return await Promise.all(
      messages.map(async (message) => {
        const id = await hashMessage(message);
        this.messages.set(id, message);
        return id;
      }),
    );
  }

  /**
   * Eviction at a step boundary (harness.md §6): keeps what the given heads reach, by reference counts
   * over segments, so a step costs a diff of segment lists, not a walk of the whole head.
   */
  keep(roots: readonly HeadRecord[]): void {
    const kept = new Set(roots.flatMap((head) => head.segments));
    for (const id of kept)
      if (!this.segmentRoots.has(id))
        for (const entry of this.segments.get(id)?.entries ?? [])
          this.entryRefs.set(entry, (this.entryRefs.get(entry) ?? 0) + 1);
    for (const id of this.segmentRoots)
      if (!kept.has(id))
        for (const entry of this.segments.get(id)?.entries ?? []) {
          const left = (this.entryRefs.get(entry) ?? 1) - 1;
          if (left > 0) this.entryRefs.set(entry, left);
          else this.entryRefs.delete(entry);
        }
    this.segmentRoots.clear();
    for (const id of kept) this.segmentRoots.add(id);
    for (const [id, segment] of this.segments)
      if (!kept.has(id)) {
        this.segments.delete(id);
        this.toSpill(`s:${id}`, segment);
      }
    const keepHeads = new Set(roots.map((head) => head.id));
    for (const [id, head] of this.heads)
      if (!keepHeads.has(id)) {
        this.heads.delete(id);
        this.toSpill(`h:${id}`, head);
      }
    for (const [id, message] of this.messages)
      if (!this.entryRefs.has(id)) {
        this.messages.delete(id);
        this.tokens.delete(id);
        if (!this.toSpill(`m:${id}`, message)) this.stored.delete(id);
      }
  }

  private unspill<T extends Message | Segment | HeadRecord>(key: string): T | undefined {
    const found = this.spill.get(key);
    if (found === undefined) return undefined;
    this.spill.delete(key);
    this.spillBytes -= found.size;
    return found.value as T;
  }

  /** Whether the record stays in the spill. */
  private toSpill(key: string, value: Message | Segment | HeadRecord): boolean {
    if (this.cacheBytes <= 0) return false;
    const size = JSON.stringify(value).length;
    this.spill.set(key, { value, size });
    this.spillBytes += size;
    for (const [oldest, { size: bytes }] of this.spill) {
      if (this.spillBytes <= this.cacheBytes) break;
      this.spill.delete(oldest);
      this.spillBytes -= bytes;
      if (oldest.startsWith('m:')) this.stored.delete(oldest.slice(2));
    }
    return this.spill.has(key);
  }

  retention(): RetentionSnapshot {
    let bytes = 0;
    for (const message of this.messages.values()) bytes += JSON.stringify(message).length;
    return { entries: this.messages.size, segments: this.segments.size, bytes, spillBytes: this.spillBytes };
  }
}

/** A context's transient head within a step: a pointer's (`scope`), or a local one (a fork, or empty). */
export interface Staged {
  readonly scope?: ScopeId;
  readonly from?: HeadRecord;
  entries: EntryId[];
  /** Hashes appended since the step started; undefined after an edit. */
  appended: EntryId[] | undefined;
}

export const local = (entries: readonly EntryId[] = []): Staged => ({
  entries: [...entries],
  appended: [],
});

/** What a top-level step will commit: calls, new messages, pointer moves, compactions, blobs. */
export class Collector {
  readonly calls: OpenCall[] = [];
  readonly messages = new Map<EntryId, Message>();
  readonly pointers = new Map<ScopeId, Staged>();
  readonly compaction = new Map<ScopeId, Compaction>();
  readonly blobs = new Map<BlobRef, Json>();
  /** Who writes each pointer in this step (harness.md §8). */
  readonly holders = new Map<ScopeId, string>();

  stage(scope: ScopeId, from: HeadRecord | undefined, entries: readonly EntryId[]): Staged {
    let staged = this.pointers.get(scope);
    if (staged === undefined) {
      staged = { scope, ...(from !== undefined ? { from } : {}), ...local(entries) };
      this.pointers.set(scope, staged);
    }
    return staged;
  }

  append(staged: Staged, ids: readonly EntryId[], messages: readonly Message[]): void {
    staged.entries.push(...ids);
    staged.appended?.push(...ids);
    ids.forEach((id, i) => {
      const message = messages[i];
      if (message !== undefined) this.messages.set(id, message);
    });
  }

  edit(staged: Staged, entries: readonly EntryId[]): void {
    staged.entries = [...entries];
    staged.appended = undefined;
  }

  /** The new messages the staged pointers name, less those known to be stored. */
  referenced(cache: WorkingSet): Message[] {
    const named = new Set<EntryId>();
    for (const staged of this.pointers.values()) for (const id of staged.appended ?? staged.entries) named.add(id);
    return [...this.messages].flatMap(([id, message]) => (named.has(id) && !cache.stored.has(id) ? [message] : []));
  }

  /** The moves, built with the shared head scheme from the working set. */
  async moves(cache: WorkingSet): Promise<Move[]> {
    const out: Move[] = [];
    for (const [pointer, staged] of this.pointers) {
      if (staged.appended !== undefined && staged.appended.length === 0) continue;
      const edit = staged.appended !== undefined ? { append: staged.appended } : { entries: staged.entries };
      const move = await buildMove({
        pointer,
        ...(staged.from !== undefined ? { from: staged.from } : {}),
        edit,
        read: cache.readSegments,
      });
      for (const segment of move.segments) cache.segments.set(segment.id, segment);
      cache.heads.set(move.head.id, move.head);
      out.push(move);
    }
    return out;
  }
}
