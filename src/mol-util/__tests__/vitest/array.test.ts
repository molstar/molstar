import { describe, it, expect } from 'vitest'; // <-- import testing functions
import { arrayMax } from '../../array';

describe('array utilities', () => {

  describe('arrayMax', () => {
    it('should find the correct element', () => {
      const arr = [10, 20, 30, 40];
      const result = arrayMax(arr);
      expect(result).toBe(40);
    });
  });
});
