# Drift: integer columns with missing values

Where: format.md and language-description.md (`column` with `missing`: "a validity list", no value shape). Lane A
refuses `missing` on integer columns (`repr.schema.language.column`).

Options: (1) a `{ values, valid }` value shape; (2) sentinel values per type. Recommended: decide when the slice first
needs one; keep the refusal until then.
