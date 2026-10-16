/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { budgetCases } from './cases/budgets.js';
import { coverageCases } from './cases/coverage.js';
import { homeCases } from './cases/homes.js';
import { journalCases } from './cases/journal.js';
import { nestedCases } from './cases/nested.js';
import { observationCases } from './cases/observation.js';
import { performanceCases } from './cases/performance.js';
import { projectionCases } from './cases/projection.js';
import { regressionCases } from './cases/regressions.js';
import { stepCases } from './cases/steps.js';
import { toolCases } from './cases/tools.js';
import { turnCases } from './cases/turns.js';
import type { ConformanceFixture } from './session.js';

/** Every conformance case, against one backend. Call it from a Jest file. */
export function runHarnessConformance(fixture: ConformanceFixture): void {
  journalCases(fixture);
  turnCases(fixture);
  homeCases(fixture);
  stepCases(fixture);
  toolCases(fixture);
  nestedCases(fixture);
  budgetCases(fixture);
  projectionCases(fixture);
  observationCases(fixture);
  performanceCases(fixture);
  coverageCases(fixture);
  regressionCases(fixture);
}
