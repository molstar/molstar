/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ReadonlyState } from '@molstar/prism-ui/core';
import type { OperationDefinition } from '../api/area.js';

// What the agent may do without asking (agent.md §4; decision 26: as little in the way as is reasonable).
// - Preauthorized: operations that change only the viewer's state while the turn can be undone. When it cannot (undo
//   is off for a large scene, sequence.md §7.4), the first changing call asks, once for the turn.
// - Network, per address: an operation says which addresses a call reaches (`addresses`); the known data services,
//   the hosts the person wrote in the request and the domains they allowed go; any other address asks, and so does an
//   operation that does not say. Only https goes without asking.
// - Ask every time: removing or replacing (`destructive`), sending out (`upload`), and a change no undo restores
//   (`undoable: false`).
// Decided when a call is made, on the values the operation will run with (the tool's `approvedBy`).

/** Network addresses the agent may reach without asking: the known data services (https only). */
export const KnownServices: readonly string[] = [
  'rcsb.org',
  'ebi.ac.uk',
  'pdbj.org',
  'alphafold.com',
  'modelarchive.org',
  'uniprot.org',
  'pdbe.org',
  // Zenodo: a public research-data archive (its records' file lists and files).
  'zenodo.org',
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

/**
 * Whether an address goes without asking: https (written whole, nothing trimmed or guessed) to a known service, a host
 * the person wrote in this turn's request, or a domain they allowed.
 */
export function allowedAddress(url: string, rules?: DomainRules, turnHosts?: ReadonlySet<string>): boolean {
  const host = hostOf(url);
  if (!host) return false;
  return KnownServices.some((d) => underDomain(host, d)) || !!turnHosts?.has(host) || !!rules?.allows(url);
}

/** The https hosts a request names (the person wrote them: the turn may reach them without asking). */
export function addressesInRequest(text: string): string[] {
  const hosts = new Set<string>();
  for (const m of text.matchAll(/https:\/\/[^\s"'<>()]+/gi)) {
    const host = hostOf(m[0]);
    if (host) hosts.add(host);
  }
  return [...hosts];
}

/** Whether a call of the operation may ask at all (a read or a view change that reaches no address never does). */
export function mayAsk(def: OperationDefinition<any, any>): boolean {
  return !!def.approval || (!def.read && (def.history !== false || def.undoable === false));
}

export interface AgentPolicy {
  /** Whether a call of this operation with these values (as it will run) goes without asking now. */
  approves(def: OperationDefinition<any, any>, values: Readonly<Record<string, unknown>>): boolean;
  /** The person allowed a call while undo was off: changing calls go for the rest of the turn. */
  allowedForTurn(): void;
  /** A turn begins: what was allowed for the last one is not for this one; the hosts its request names are. */
  beginTurn(requestHosts?: readonly string[]): void;
  readonly rules: DomainRules;
  /**
   * The https hosts a network call reaches that are not allowed yet (what "Allow this domain" would allow); empty when
   * it reaches none, or asks for another reason (an address not https, an operation that does not say).
   */
  blockedHosts(def: OperationDefinition<any, any>, values: Readonly<Record<string, unknown>>): string[];
  /** A tool's blocked hosts from its arguments (registered by the tools that reach addresses). */
  registerHosts(tool: string, hosts: (args: unknown) => string[]): void;
  /** A waiting call's blocked hosts, by its tool and arguments; empty when the tool reaches none. */
  hostsOf(tool: string, args: unknown): string[];
}

export function createAgentPolicy(o: { checkpointing: ReadonlyState<boolean>; rules?: DomainRules }): AgentPolicy {
  const rules = o.rules ?? createDomainRules();
  let turnAllowed = false;
  let turnHosts = new Set<string>();
  const tools = new Map<string, (args: unknown) => string[]>();
  return {
    rules,
    blockedHosts(def, values) {
      if (def.approval !== 'network' || !def.addresses) return [];
      const hosts = new Set<string>();
      for (const a of def.addresses(values as never)) {
        if (typeof a !== 'string' || allowedAddress(a, rules, turnHosts)) continue;
        const host = hostOf(a);
        // Not https (or not a URL written whole): no domain rule allows it; it asks every time.
        if (!host) return [];
        hosts.add(host);
      }
      return [...hosts];
    },
    registerHosts(tool, hosts) {
      tools.set(tool, hosts);
    },
    hostsOf(tool, args) {
      try {
        return tools.get(tool)?.(args) ?? [];
      } catch {
        return [];
      }
    },
    approves(def, values) {
      if (def.approval === 'destructive' || def.approval === 'upload') return false;
      if (def.approval === 'network') {
        // An operation that does not say what it reaches always asks.
        if (!def.addresses) return false;
        const addresses = def.addresses(values as never);
        if (!addresses.every((a) => typeof a === 'string' && allowedAddress(a, rules, turnHosts))) return false;
      }
      if (def.read) return true;
      // What no undo restores is the person's to allow (agent.md §4: preauthorized only while undoable): a snapshot
      // added, a file made for them (an export, a saved session), whether or not it is a history step.
      if (def.undoable === false) return false;
      if (def.history === false) return true;
      // The turn cannot be undone (undo off for a large scene): the first change asks, once for the turn.
      if (!o.checkpointing.value && !turnAllowed) return false;
      return true;
    },
    allowedForTurn: () => {
      // Only the person's allow given while the turn cannot be undone waives that ask (an allowed fetch does not).
      if (!o.checkpointing.value) turnAllowed = true;
    },
    beginTurn: (requestHosts = []) => {
      turnAllowed = false;
      turnHosts = new Set(requestHosts);
    },
  };
}
