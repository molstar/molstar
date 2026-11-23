# Records

One record per Seq\* execution session, kept permanently (§permanent-record). A record says **how we got here**: what
the author asked, what was decided, what was delivered and how it was checked, what was learned, and what was left
open. The designs in `.seq/design/` say **what is true now**; a record never replaces them.

## Layout

```
.seq/records/
  <work-id>/
    <yyyy-mm-dd>-<short-name>.md   # one file per execution session, named by the day it ended
```

The work id is the one used in `.seq/execution/<work-id>/`. A record folder holds the record file only; the work folder
stays in `.seq/execution/` after the close, never moved or deleted unless the author says so.

## Rules

- **Written at the end of a session** (and at a compaction or handoff), from the work folder. Never edited afterwards;
  a correction is a new record that says so.
- **Decisions say how they were made** (§authors-words-decide): one the author spelled out, with the quote; one
  accepted as recommended, marked so, with the alternatives not taken and the reason.
- **Short**: one to two pages. Commits by hash.
- **Decisions are copied, not moved**: each goes into the design it belongs to as well.
- **Machine-invariant** (§machine-invariant): no local paths; other repositories named in words.

## Template

The session record in [workflow-files](../oz/_core/patterns/workflow-files.md): Asked, Decided, Delivered, Learned,
Left open, and Self-check when the session changed the principles.
