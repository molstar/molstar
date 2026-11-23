# Lane EO: the third slice's edge-case outcomes, drafted (read-only)

**Goal**: SR found `.scaffolds/seq/design/edge-cases.md` has no "Outcomes of the third slice" (`reconcile.md`, list 4). Draft it,
so finalize pastes a checked table. **Read only**: your report is the table.

For each EC cited in `packages/seq` with no outcome row (SR's list: 11, 20, 33, 35, 38–41, 43, 44, 50, 51, 57, 59,
62–65, 67, 87, 92, 94–97, 99–103, 106–109, 112, 113, 115–118, 128, 129, 132, 137, 138, 145, 195, 204, 206, 209, 232,
233; and EC-85, EC-126, EC-130, EC-135 to recheck): read its text in `edge-cases.md`, find the test that shows the
behavior (by an `EC-<n>` citation, else by behavior: grep the tests for what the EC describes), run that one test
file (`pnpm test:unit <file>`), and write the row in the existing tables' shape: `| EC | outcome (kept, changed:
how, not met: what is missing) | the test (package-relative path and test name) |`. An EC whose behavior has no test:
"not met" with what a test would need, never a test that does not exist. Mark EC-154/155 (after RK) and EC-169–175
(after CS) "pending".

**Model**: fast (sonnet). **Where**: the main checkout, read only: no file written, no installs or builds, no git
writes; `pnpm test:unit <one file>` is fine.

**Stop when** the table is done, or after about an hour: report what you have.

**Report**: the table, ready to paste under "## Outcomes of the third slice"; the count kept / changed / not met;
any EC whose text the code contradicts (a design question, not a row).
