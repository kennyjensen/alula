// SPDX-License-Identifier: GPL-2.0-or-later
// Detached diagnostic data only. Reporting a failed solve must neither lose
// a typed inner error nor fail on cycles, getters, or non-JSON values.
export function serializeDiagnosticError(error) {
  const seenCauses = new Set(), stack = new Set();
  let remaining = 20000;
  const read = (object, key) => {
    try { return object[key]; } catch { return '[Unreadable property]'; }
  };
  const plain = (value, depth = 0) => {
    if (--remaining < 0) return '[Truncated: diagnostic budget]';
    if (value === null || typeof value === 'boolean') return value;
    if (value === undefined) return undefined;
    if (typeof value === 'number') return Number.isFinite(value) ? value : String(value);
    if (typeof value === 'string') return value.length <= 32768 ? value : `${value.slice(0, 32768)}[Truncated: string]`;
    if (typeof value === 'bigint' || typeof value === 'symbol') return String(value);
    if (typeof value === 'function') return '[Function]';
    if (stack.has(value)) return '[Circular]';
    if (depth >= 16) return '[Truncated: diagnostic depth]';
    stack.add(value);
    try {
      if (Array.isArray(value) || ArrayBuffer.isView(value)) {
        if (value instanceof DataView) return plain(new Uint8Array(value.buffer, value.byteOffset, value.byteLength), depth + 1);
        const length = read(value, 'length'), count = Math.min(length, 1024);
        if (!Number.isInteger(length) || length < 0) return '[Unreadable array]';
        const result = Array.from({ length: count }, (_, i) => plain(read(value, i), depth + 1) ?? null);
        if (count < length) result.push(`[Truncated: ${length - count} array entries]`);
        return result;
      }
      let keys;
      try { keys = Object.keys(value); } catch { return '[Unreadable object]'; }
      const entries = keys.slice(0, 256).flatMap(key => {
        const copied = plain(read(value, key), depth + 1);
        return copied === undefined ? [] : [[key, copied]];
      });
      if (keys.length > 256) entries.push(['diagnosticTruncation', `${keys.length - 256} object entries omitted`]);
      return Object.fromEntries(entries);
    } catch {
      return '[Unreadable diagnostic value]';
    } finally { stack.delete(value); }
  };
  const cause = (value, depth = 0) => {
    if (value === null || typeof value !== 'object') return plain(value);
    if (seenCauses.has(value)) return '[Circular cause]';
    if (depth >= 8) return '[Truncated: cause depth]';
    seenCauses.add(value);
    const result = Object.fromEntries(['name', 'message', 'code', 'stage', 'reason', 'diagnostics', 'initialization', 'stack']
      .flatMap(key => {
        const copied = plain(read(value, key));
        return copied === undefined ? [] : [[key, copied]];
      }));
    const nested = read(value, 'cause');
    if (nested !== undefined) result.cause = cause(nested, depth + 1);
    return result;
  };
  const result = cause(error);
  return result && typeof result === 'object' ? result : { message: result ?? 'Unknown error' };
}
