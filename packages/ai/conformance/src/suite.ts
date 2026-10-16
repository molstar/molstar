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

/** Every conformance case, on the in-memory harness (testing.md). Call it from a test file. */
export function runHarnessConformance(): void {
  journalCases();
  turnCases();
  homeCases();
  stepCases();
  toolCases();
  nestedCases();
  budgetCases();
  projectionCases();
  observationCases();
  performanceCases();
  coverageCases();
  regressionCases();
}
