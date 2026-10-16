/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Message, ModelConfiguration } from '../../messages/index.js';
import { type DriveContext, type HeadEdit, MAIN } from '../../harness/index.js';
import type { CompactionDomain } from './domain.js';
import { type CompactionPlan, editTokens, type PlanInput, planCompaction, sizeOf } from './plan.js';
import { carryOverInstruction, type CompactionSettings } from './settings.js';

const textMessage = (role: 'user' | 'system', text: string): Message => ({ role, parts: [{ kind: 'text', text }] });

export interface CompactionOptions {
  settings: CompactionSettings;
  domain: CompactionDomain;
  /** The model that writes the carry-over (the turn's, or at a turn's start the coming turn's). */
  model: ModelConfiguration;
  /** force: the context no longer fits; trigger: mid-turn, at the trigger; soft: a turn's start. */
  when: 'force' | 'trigger' | 'soft';
  /** The context's size now (the decision's), when it has one. */
  size?: number;
  /** The size the last compaction left, if any: no compaction again until it has grown by the headroom. */
  lastSize?: number;
  onPlan?: (plan: CompactionPlan) => void;
}

/**
 * Whether the context needs condensing now, and the new head if so (session.md §4): one nested model call
 * on the context as it stands asks what to carry over; the plan builds the new head from its answer (or from a
 * digest, when the call fails or the context cannot be sent any more). Above its target the plan is accepted all the
 * same, with a size note for the next one.
 */
export async function runCompaction(context: DriveContext, options: CompactionOptions): Promise<HeadEdit | undefined> {
  if (context.depth > 0) return undefined;
  const { settings, when, domain } = options;
  const trigger = settings.triggerInputTokens;
  // After a compaction that stayed above its target, wait for the usual headroom before the next.
  const headroom = trigger - settings.targetContextTokens;
  const threshold = Math.max(trigger, options.lastSize !== undefined ? options.lastSize + headroom : 0);
  const sizeNow = options.size ?? 0;
  if (when === 'trigger' && sizeNow < threshold) return undefined;
  if (when === 'soft' && sizeNow < settings.softShare * threshold) return undefined;
  const { messages } = await context.journal.read(MAIN);
  const entries = context.head.entries.map((id, i) => ({ id, message: messages[i] as Message }));
  if (when === 'soft' && !entries.some((e) => domain.isRequest(e.message))) return undefined;
  const midTurn = when === 'trigger';
  let answer: string | undefined;
  let call: CompactionPlan['report']['call'];
  let callFailure: string | undefined;
  if (when !== 'force') {
    const model = options.model;
    const result = await context
      .run({
        input: [...messages, textMessage('user', carryOverInstruction(midTurn))],
        tools: { enabled: [], preauthorized: [], definitions: [] },
        model: { ...model, parameters: { ...model.parameters, maxOutputTokens: settings.maxOutputTokens } },
        limits: { modelCalls: 1 },
      })
      .catch((e: unknown) => {
        callFailure = e instanceof Error ? e.message : String(e);
        return undefined;
      });
    if (result !== undefined && result.status !== 'completed')
      callFailure = JSON.stringify(result.failure ?? result.status);
    if (result?.status === 'completed' && result.output !== undefined) {
      const parts = result.output.parts;
      if (!parts.some((p) => p.kind === 'tool-call'))
        answer = parts.map((p) => (p.kind === 'text' ? p.text : '')).join('');
      const u = result.usage;
      call = {
        input: u.uncachedInput + u.cacheRead + u.cacheWrite,
        cached: u.cacheRead,
        output: u.output,
        reasoning: u.reasoning,
      };
    }
  }
  const target = settings.targetContextTokens;
  const pins = domain.pins?.() ?? { skills: [], tools: [] };
  const memory = domain.memory?.();
  const make = (lean: boolean, sizeNote?: PlanInput['sizeNote']) =>
    planCompaction({
      entries,
      domain,
      current: when !== 'soft',
      answer,
      pins: lean ? { skills: [], tools: pins.tools } : pins,
      settings: lean ? { ...settings, recentTurnsTokens: 0, maxSummaryTokens: 0 } : settings,
      ...(memory ? { memory } : {}),
      ...(sizeNote ? { sizeNote } : {}),
    });
  const oldTokens = options.size ?? entries.reduce((a, e) => a + sizeOf(e.message), 0);
  const measure = (p: CompactionPlan) => editTokens(p.edit, entries, oldTokens);
  // Aim at the target: try a leaner plan (no recent turns, no skills again); still above, ask the next compaction
  // to condense harder, but only when the carried-over part is what is over.
  let lean = false;
  let plan = make(false);
  let tokens = measure(plan);
  if (tokens > target) {
    const leaner = make(true);
    const t = measure(leaner);
    if (t < tokens) [plan, tokens, lean] = [leaner, t, true];
  }
  let fixed: number | undefined;
  if (tokens > target) {
    const leanest = lean ? plan : make(true);
    const notes = new Set(leanest.edit.entries.filter((e) => typeof e !== 'string' && e.role === 'user'));
    fixed = editTokens({ entries: leanest.edit.entries.filter((e) => !notes.has(e as never)) }, entries, oldTokens);
    if (fixed < target) {
      plan = make(lean, { estimated: tokens, target, fixed });
      tokens = measure(plan);
    }
  }
  plan.report.estimatedTokens = tokens;
  if (tokens > target) plan.report.overTarget = true;
  if (fixed !== undefined) {
    plan.report.fixedTokens = fixed;
    if (fixed >= target) plan.report.targetOutOfReach = true;
  }
  if (call) plan.report.call = call;
  if (callFailure !== undefined) plan.report.callFailure = callFailure;
  options.onPlan?.(plan);
  return plan.edit;
}
