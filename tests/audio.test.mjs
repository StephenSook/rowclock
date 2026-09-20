// Cross-instrument arithmetic. These are closed-form checks: the acoustic
// reference and the implied row time are related by one equation, and the
// plausibility window is a physical fact about consumer sensors, not a
// tolerance chosen to make a result pass.

import { test } from 'node:test';
import assert from 'node:assert/strict';
import { lineTimeFromAcousticReference, crossCheck } from '../src/audio.js';

test('the row time falls out of the two instruments', () => {
  // A drill turning at 333 rev/s with two marks sounds at 667 Hz.
  // The camera, reading a row every 25.5 us, sees that as:
  const trueLineTime = 25.5e-6;
  const acoustic = 667;
  const cyclesPerRow = acoustic * trueLineTime;

  const solved = lineTimeFromAcousticReference(cyclesPerRow, acoustic);
  assert.ok(Math.abs(solved - trueLineTime) < 1e-12, `solved ${solved}`);
  assert.ok(Math.abs(solved * 1e6 - 25.5) < 1e-6);
});

test('it refuses impossible inputs rather than returning a number', () => {
  assert.throws(() => lineTimeFromAcousticReference(0.017, 0), /positive/);
  assert.throws(() => lineTimeFromAcousticReference(0, 667), /positive/);
});

test('crossCheck accepts a physically plausible pairing', () => {
  const c = crossCheck(667 * 25.5e-6, 667);
  assert.equal(c.agree, true);
  assert.ok(Math.abs(c.lineTimeUs - 25.5) < 0.01, `row time ${c.lineTimeUs}`);
  assert.ok(Math.abs(c.rowRateHz - 39215) < 50, `row rate ${c.rowRateHz}`);
});

test('crossCheck REJECTS two numbers that agree into an impossible sensor', () => {
  // The dangerous case: the camera locked onto a harmonic, so the ratio is
  // tidy and the answer is nonsense. A row time of 255 us would mean a
  // 1080-row frame takes 275 ms to read out, which is 8 frames per second.
  const bad = crossCheck(667 * 255e-6, 667);
  assert.equal(bad.agree, false);
  assert.match(bad.reason, /outside the/);
});

test('crossCheck rejects a row time that is far too fast to be consumer hardware', () => {
  const bad = crossCheck(667 * 0.5e-6, 667);
  assert.equal(bad.agree, false);
});

test('crossCheck offers the harmonic readings that are also plausible', () => {
  // A five-blade fan at 12 rev/s sounds at 60 Hz, and the camera may lock to
  // the shaft rate instead. Surfacing both is honest; silently picking one is
  // not.
  const c = crossCheck(60 * 30e-6, 60);
  assert.equal(c.agree, true);
  assert.ok(c.harmonicCandidates.length >= 1);
  for (const h of c.harmonicCandidates) {
    assert.ok(h.us >= 4 && h.us <= 80, `implausible candidate offered: ${h.us}`);
  }
});

test('crossCheck says so when an instrument has no reading', () => {
  assert.equal(crossCheck(0, 667).agree, false);
  assert.equal(crossCheck(0.017, 0).agree, false);
  assert.match(crossCheck(0, 667).reason, /no reading/);
});
