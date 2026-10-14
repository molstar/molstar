/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The adapter's issues, each at its code's registered severity (seqspec:format#validation): Nightingale reports
// kernel codes only, so the kernel's code table decides.
import { createRegistries, issueOf, KernelTables, register, type SvsRegistries } from "@molstar/seq-core/format";
import type { Issue } from "@molstar/seq-core/model";

let codes: SvsRegistries | null = null;

/** An issue under a kernel code, its severity the registered one's. */
export function nightingaleIssue(code: string, path: string, message: string): Issue {
  if (!codes) {
    codes = createRegistries();
    register(codes, { issueCodes: KernelTables.issueCodes });
  }
  return issueOf(codes, code, path, message);
}
