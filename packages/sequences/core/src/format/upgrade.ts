/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `upgrade` (seqspec:format#versioning): one step per older version, up to the version this package reads.

import { SvsFormat, SvsVersion } from "./doc.ts";
import { isPlainObject } from "./schema/field-schema.ts";

/** A document of version `from` (`major.minor`) to `to`; `run` gets a copy it may change. */
export interface UpgradeStep {
  readonly from: string;
  readonly to: string;
  run(doc: Record<string, unknown>): Record<string, unknown>;
}

/** The step registry: none yet, 1.0 being the first version. */
export const UpgradeSteps: readonly UpgradeStep[] = [];

const minorOf = (v: unknown) => (typeof v === "string" ? (/^(\d+\.\d+)\.\d+$/.exec(v)?.[1] ?? null) : null);

/**
 * Runs the step from the document's version while there is one and it is older than `SvsVersion`; each step sets
 * the version to its `to`. Anything that is not a SeqViewSpec, or needs no step, is returned as it is.
 */
export function upgrade(doc: unknown, steps: readonly UpgradeStep[] = UpgradeSteps): { doc: unknown; steps: string[] } {
  const done: string[] = [];
  if (!isPlainObject(doc) || doc.format !== SvsFormat) return { doc, steps: done };
  let d = doc;
  for (let v = minorOf(d.version); v !== minorOf(SvsVersion) && done.length < steps.length; v = minorOf(d.version)) {
    const step = steps.find((s) => s.from === v);
    if (!step) break;
    d = { ...step.run(structuredClone(d)), version: `${step.to}.0` };
    done.push(`${step.from} -> ${step.to}`);
  }
  return { doc: d, steps: done };
}
