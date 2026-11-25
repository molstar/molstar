/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CommandFormat, FormatCheck } from './types.js';

/**
 * A format's check, guarded (commands §2.2): one that throws reads as text the format does not take, with what it threw
 * as the message. Why: a format is the host's code, and every reader (typing, the tail, Enter, the run's last check)
 * must stay in a state the person can edit; a throw out of the run's check left the invocation running for good.
 */
export function checkFormat<Ctx>(format: CommandFormat<Ctx>, text: string, ctx: Ctx): FormatCheck {
  try {
    return format.check(text, ctx);
  } catch (e) {
    return { ok: false, message: e instanceof Error ? e.message : String(e) };
  }
}
