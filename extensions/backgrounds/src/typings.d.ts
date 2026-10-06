/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

declare module '*.jpg' {
  const value: string;
  // biome-ignore lint/style/noDefaultExport: Image loaders expose the asset URL as a default export.
  export { value as default };
}
