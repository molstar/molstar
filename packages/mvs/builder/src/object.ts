/** Small object helpers used by the standalone MVS builder. */
export function isPlainObject(value: any): boolean {
  return typeof value === 'object' && value !== null && !Array.isArray(value);
}

export function mapObjectMap<T, S>(obj: { [key: string]: T }, f: (value: T) => S): { [key: string]: S } {
  const result: { [key: string]: S } = {};
  for (const key of Object.keys(obj)) result[key] = f(obj[key]);
  return result;
}

export function omitObjectKeys<T extends {}, K extends keyof T>(obj: T, keys: readonly K[]): Omit<T, K> {
  const result: T = { ...obj };
  for (const key of keys) delete result[key];
  return result as Omit<T, K>;
}

export function pickObjectKeys<T extends {}, K extends keyof T>(obj: T, keys: readonly K[]): Pick<T, K> {
  const result: Partial<Pick<T, K>> = {};
  for (const key of keys) if (Object.hasOwn(obj, key)) result[key] = obj[key];
  return result as Pick<T, K>;
}

export function deepClone<T>(source: T): T {
  if (source === null || typeof source !== 'object') return source;
  if (Array.isArray(source)) return source.map(deepClone) as T;
  if (!('prototype' in source)) {
    const copy: Record<string, unknown> = {};
    for (const key in source)
      if (Object.prototype.hasOwnProperty.call(source, key)) copy[key] = deepClone((source as any)[key]);
    return copy as T;
  }
  throw new Error(`Can't clone, type "${typeof source}" unsupported`);
}
