// SPDX-License-Identifier: GPL-2.0-or-later
// Numerical termination is distinct from a user observer or arithmetic error.
class SlorNumericalError extends Error {
  constructor(message, termination) { super(message); this.termination = termination; }
}
export const slorNumericalError = (message, termination) => new SlorNumericalError(message, termination);
export const slorErrorTermination = error => error instanceof SlorNumericalError
  ? { origin: 'solver', termination: error.termination }
  : { origin: error?.code === 'slor-observer-failed' ? 'observer' : 'exception', termination: 'exception' };
export function slorObserver(observer, history, nodes) {
  if (observer === undefined) return;
  // The smoother has already detached this history row and grid before
  // invoking the initializer. Its recorder runs before this ownership transfer
  // and copies any state it retains. A second deep copy only adds allocation.
  try { observer(history, nodes); }
  catch (cause) { throw Object.assign(new Error(cause?.message ?? String(cause), { cause }), { code: 'slor-observer-failed' }); }
}
export function tagSlorTermination(result, termination) {
  // Keep successful legacy reports and nested diagnostic shapes unchanged.
  // The initializer explicitly serializes this provenance when it invokes
  // partial-guess recovery; it is never inferred from an error's wording.
  if (!result.converged) Object.defineProperty(result, 'termination', { value: termination });
  return result;
}
