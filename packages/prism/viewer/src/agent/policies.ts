/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json } from '@molstar/ai-agent';
import type { ReadonlyState } from '@molstar/prism-ui/core';
import type { OperationDefinition } from '../api/area.js';

// What the agent may do without asking (agent.md §4; decision 26: as little in the way as is reasonable).
// - Preauthorized: operations that change only the viewer's state while the turn can be undone. When it cannot (undo
//   is off for a large scene, sequence.md §7.4), the first changing call asks, once for the turn.
// - Network, per address: the known data services and domains the person allowed go; any other address asks.
// - Ask every time: removing or replacing (`destructive`) and sending out (`upload`).
// Decided when a call is made (the session's permission decision runs the tool's `approvedBy`).

/** Network addresses the agent may reach without asking: the known data services (https only). */
export const KnownServices: readonly string[] = [
  'rcsb.org',
  'ebi.ac.uk',
  'pdbj.org',
  'alphafold.com',
  'modelarchive.org',
  'uniprot.org',
  'pdbe.org',
];

/** The host of an https address, or undefined. */
export function hostOf(url: string): string | undefined {
  try {
    const u = new URL(url);
    return u.protocol === 'https:' ? u.hostname.toLowerCase() : undefined;
  } catch {
    return undefined;
  }
}

const underDomain = (host: string, domain: string) => host === domain || host.endsWith(`.${domain}`);

/** Domains the person allowed (Allow this domain): for the session, or kept by the host (in this browser). */
export interface DomainRules {
  allows(url: string): boolean;
  allow(domain: string, scope: 'session' | 'browser'): void;
  remove(domain: string): void;
  list(): readonly { readonly domain: string; readonly scope: 'session' | 'browser' }[];
}

/** Domain rules in memory, those kept `browser` handed to `save` (the assistant keeps them; agent.md §4). */
export function createDomainRules(
  kept: readonly string[] = [],
  save?: (browser: readonly string[]) => void,
): DomainRules {
  const rules = new Map<string, 'session' | 'browser'>(kept.map((d) => [d, 'browser']));
  const persist = () => save?.([...rules].filter(([, s]) => s === 'browser').map(([d]) => d));
  return {
    allows: (url) => {
      const host = hostOf(url);
      return !!host && [...rules.keys()].some((d) => underDomain(host, d));
    },
    allow(domain, scope) {
      rules.set(domain.toLowerCase(), scope);
      if (scope === 'browser') persist();
    },
    remove(domain) {
      if (rules.delete(domain.toLowerCase())) persist();
    },
    list: () => [...rules].map(([domain, scope]) => ({ domain, scope })),
  };
}

/** Whether an address goes without asking: https to a known service, or a domain the person allowed. */
export function allowedAddress(url: string, rules?: DomainRules): boolean {
  const host = hostOf(url);
  if (!host) return false;
  return KnownServices.some((d) => underDomain(host, d)) || !!rules?.allows(url);
}

/** The addresses a call's arguments name (a field whose value is one). */
export function addressesIn(args: Json): string[] {
  if (!args || typeof args !== 'object' || Array.isArray(args)) return [];
  return Object.values(args).filter((v): v is string => typeof v === 'string' && /^[a-z]+:\/\//i.test(v));
}

export interface AgentPolicy {
  /** Whether a call of this operation with these arguments goes without asking now. */
  approves(def: OperationDefinition<any, any>, args: Json): boolean;
  /** The person allowed a call this turn: changing calls go for the rest of it, undo off or not. */
  allowedForTurn(): void;
  /** A turn begins: what was allowed for the last one is not for this one. */
  beginTurn(): void;
  readonly rules: DomainRules;
}

export function createAgentPolicy(o: { checkpointing: ReadonlyState<boolean>; rules?: DomainRules }): AgentPolicy {
  const rules = o.rules ?? createDomainRules();
  let turnAllowed = false;
  return {
    rules,
    approves(def, args) {
      if (def.approval === 'destructive' || def.approval === 'upload') return false;
      if (def.approval === 'network' && !addressesIn(args).every((a) => allowedAddress(a, rules))) return false;
      const changes = !def.read && def.history !== false && def.undoable !== false;
      // The turn cannot be undone (undo off for a large scene): the first change asks, once for the turn.
      if (changes && !o.checkpointing.value && !turnAllowed) return false;
      return true;
    },
    allowedForTurn: () => {
      turnAllowed = true;
    },
    beginTurn: () => {
      turnAllowed = false;
    },
  };
}
