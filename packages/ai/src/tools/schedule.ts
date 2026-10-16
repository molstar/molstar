/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * The order a model step's tool calls run in (the drive's schedule decision, harness.md §7): stages run one after
 * another, the calls in a stage together. A run of consecutive reads (read-only tools: their result depends on their
 * arguments and the state alone) is one stage; every other call (an edit, a view command, a disclosure) is a stage of
 * its own, in the model's order. So edits apply in the order the model wrote them, a read after an edit sees it, and
 * reads between edits do not wait for each other.
 */
export function toolSchedule(
  calls: readonly { callId: string; tool: { id: string }; order?: number }[],
  reads: (toolId: string) => boolean,
): string[][] {
  const ordered = [...calls].sort((a, b) => (a.order ?? 0) - (b.order ?? 0));
  const stages: string[][] = [];
  let together: string[] | undefined;
  for (const c of ordered) {
    if (reads(c.tool.id)) {
      if (!together) {
        together = [];
        stages.push(together);
      }
      together.push(c.callId);
    } else {
      together = undefined;
      stages.push([c.callId]);
    }
  }
  return stages;
}
