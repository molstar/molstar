// The defaults-designed check (the rule "defaults are designed with a picture", from Seq*'s fold-backs): a design
// under `<instance>/design/` whose `## Consequences` section names something drawn (the instance's words) has a
// `## Defaults` heading (what a reader sees with no setting). `--worklist` prints every design that fails, exempt or not, so the
// close can see what is owed.
// The method's, not an instance's: any scaffold whose designs add something drawn owes its defaults; one that draws
// nothing does not run it. The instance's part is `<instance>/_checks/defaults-designed.json`, required:
// `{ "drawn": "<regex: the words of a drawn thing in Consequences>", "exempt": ["<design name without .md>", ...] }`;
// `exempt` lists the designs that failed when the check was written and only shrinks, a design leaves it when it gets
// its Defaults section (the check says so).
// Usage: defaults-designed.ts [root] --instance .scaffolds/<name> [--worklist]
import { existsSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { finish, instanceFile, methodArgs, section, trackedFiles, usage } from "./lib.ts";

const { root: ROOT, instance: INSTANCE, flags } = methodArgs({ flags: ["--worklist"] });
const WORKLIST = flags.has("--worklist");
const DESIGNS = `${INSTANCE}/design/`;
const { drawn, exempt = [] } = instanceFile<{ drawn?: string; exempt?: string[] }>(
  ROOT,
  INSTANCE,
  "defaults-designed.json",
);
if (!drawn) usage(`${INSTANCE}/_checks/defaults-designed.json: "drawn" (a regex) is required`);
/** Designs without a Defaults heading when the check was written. Only ever shrinks. */
const EXEMPT = new Set(exempt);

const NAMES_A_DRAWN_THING = new RegExp(drawn, "i");
const problems: string[] = [];
const worklist: string[] = [];
for (const f of trackedFiles(ROOT)) {
  if (!f.startsWith(DESIGNS) || !f.endsWith(".md") || !existsSync(join(ROOT, f))) continue;
  const id = f.slice(DESIGNS.length, -3);
  const text = readFileSync(join(ROOT, f), "utf8");
  const consequences = section(text, /^Consequences\b/);
  if (consequences === undefined || !NAMES_A_DRAWN_THING.test(consequences)) continue;
  if (section(text, /^Defaults\b/) !== undefined) continue;
  worklist.push(f);
  if (!EXEMPT.has(id))
    problems.push(`${f}: its Consequences name something drawn (${drawn}) but it has no "## Defaults" heading`);
}
if (WORKLIST)
  console.log(`defaults-designed: designs with no Defaults heading\n${worklist.map((w) => `  - ${w}`).join("\n")}`);
for (const id of EXEMPT)
  if (!worklist.includes(`${DESIGNS}${id}.md`)) problems.push(`${DESIGNS}${id}.md: passes now: drop it from defaults-designed.json`);
finish("defaults-designed", problems);
