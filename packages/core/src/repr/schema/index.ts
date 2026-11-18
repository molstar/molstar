/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/core/repr/schema: languages as data (kinds, parents, typed params including columns and tables, defaults,
// docs), validation with every issue, normalization on read, and what is generated from a description (_spec/README.md).

export {
  type CheckLanguageOptions,
  checkLanguage,
  defineLanguage,
  type LanguageCode,
  LanguageIssues,
} from './check.js';
export type { Changes, ParamChange, Where } from './changes.js';
export { languageToMarkdown, signature, withDefaults } from './docs.js';
export {
  type ArgSpec,
  type ColumnSpec,
  type ColumnValue,
  type DataSpec,
  dataFormsOf,
  type Fields,
  type KindSpec,
  type Language,
  type ParamSpec,
  type ParamsOf,
  type ParamType,
  type ParamValue,
  type ParamValues,
  type RefType,
  type TableSpec,
  type UnitsHook,
} from './language.js';
export {
  type ColumnPlaces,
  declaredColumnPlaces,
  type NormalizeOptions,
  normalize,
  normalizeValue,
  validateAndNormalize,
} from './normalize.js';
export {
  refHint,
  SchemaIssues,
  typeName,
  unionNodeMember,
  type ValidateOptions,
  validate,
  validateChanges,
  validateValue,
} from './validate.js';
