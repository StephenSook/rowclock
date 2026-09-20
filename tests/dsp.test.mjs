// Run with:  node --test tests/
//
// These tests check the solver against closed-form answers we can write down
// in advance. The synthetic signals here are a UNIT TEST of an algorithm, not
// evidence about the product: nothing in this file is ever quoted as a
// measurement of a real scene. Every product claim is measured on a real
// camera frame and is reported separately.

import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  fft,
  nextPow2,
  detrend,
  spectrum,
  dominantPeak,
  localBaseline,
  spectrumComplex,
  coherence,
  coherentPeak,
  cyclesInRecord,
  prominentPeak,
  detrendWindowFor,
  PeakTracker,
  analyseProfile,
  lineTimeFromKnownSource,
  hzFromCyclesPerRow,
  aliasNote,
} from '../src/dsp.js';

/** A sinusoid at a known cycles-per-sample rate, optionally over a scene ramp. */
function synth(n, cyclesPerSample, { amp = 1, phase = 0, dc = 0, ramp = 0, noise = 0, seed = 1 } = {}) {
  // Deterministic LCG so a failure is reproducible.
  let s = seed >>> 0;
  const rnd = () => ((s = (s * 1664525 + 1013904223) >>> 0) / 4294967296) - 0.5;
  const x = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    x[i] = dc
      + ramp * (i / n)
      + amp * Math.sin(2 * Math.PI * cyclesPerSample * i + phase)
      + noise * rnd();
  }
  return x;
}

test('nextPow2 rounds up to a power of two', () => {
  assert.equal(nextPow2(1), 1);
  assert.equal(nextPow2(1080), 2048);
  assert.equal(nextPow2(2048), 2048);
  assert.equal(nextPow2(2049), 4096);
});

test('fft rejects a non-power-of-two length rather than returning nonsense', () => {
  assert.throws(() => fft(new Float64Array(6), new Float64Array(6)), /power of two/);
});

test('fft of a unit impulse is flat at magnitude 1', () => {
  const n = 64;
  const re = new Float64Array(n);
  const im = new Float64Array(n);
  re[0] = 1;
  fft(re, im);
  for (let k = 0; k < n; k++) {
    assert.ok(Math.abs(Math.hypot(re[k], im[k]) - 1) < 1e-12, `bin ${k} was not 1`);
  }
});

test('fft puts a pure tone in exactly the bin it belongs in', () => {
  const n = 256;
  const k0 = 17;
  const re = new Float64Array(n);
  const im = new Float64Array(n);
  for (let i = 0; i < n; i++) re[i] = Math.cos((2 * Math.PI * k0 * i) / n);
  fft(re, im);
  let best = 1;
  for (let k = 1; k < n / 2; k++) if (Math.hypot(re[k], im[k]) > Math.hypot(re[best], im[best])) best = k;
  assert.equal(best, k0);
});

test('detrend removes DC and a scene ramp but keeps the tone', () => {
  const n = 1080;
  const cps = 0.05;
  const clean = synth(n, cps);
  const dirty = synth(n, cps, { dc: 120, ramp: 60 });

  const dt = detrend(dirty, 65);

  // A centered moving-average high-pass cannot drive the mean to exactly zero,
  // because the window is truncated at the two edges of the array and is
  // therefore asymmetric there. For a ramp of height `ramp` the edge residual
  // is of order ramp*half/n, so asserting 1e-6 asserts a perfection this
  // filter structurally cannot deliver. The property that actually matters is
  // that the residual is negligible NEXT TO THE TONE we are trying to keep.
  const mean = dt.reduce((a, b) => a + b, 0) / n;
  const toneAmplitude = 1; // synth default
  assert.ok(
    Math.abs(mean) < 0.01 * toneAmplitude,
    `residual DC ${mean} is not small relative to the tone`,
  );

  // The property a high-pass actually guarantees is spectral, not statistical:
  // the scene's energy lives in the lowest bins, and after detrending those
  // bins must be small next to the tone's own bin. Comparing means is the
  // wrong instrument here, because the interior legitimately still contains
  // the tone we set out to keep, so its mean is not a measure of cleanliness.
  const before = spectrum(Float64Array.from(dirty));
  const after = spectrum(dt);
  const toneBin = Math.round(cps * after.nfft);

  const lowBandEnergy = (s) => {
    let e = 0;
    for (let k = 1; k < 4; k++) e += s.mag[k] * s.mag[k];
    return e;
  };

  const ratioBefore = lowBandEnergy(before) / (before.mag[toneBin] ** 2);
  const ratioAfter = lowBandEnergy(after) / (after.mag[toneBin] ** 2);
  assert.ok(
    ratioAfter < ratioBefore / 100,
    `detrend only cut low-band energy from ${ratioBefore.toExponential(2)} to ${ratioAfter.toExponential(2)}`,
  );

  // The tone should survive with most of its amplitude intact.
  const pClean = analyseProfile(clean);
  const pDirty = analyseProfile(dirty);
  assert.ok(pClean && pDirty);
  assert.ok(
    Math.abs(pClean.cyclesPerSample - pDirty.cyclesPerSample) < 1e-4,
    `scene content moved the estimate: ${pClean.cyclesPerSample} vs ${pDirty.cyclesPerSample}`,
  );
});

test('dominantPeak recovers a known frequency to better than a tenth of a bin', () => {
  // Deliberately off-bin so the parabolic interpolation is doing real work.
  const n = 1080;
  const nfft = nextPow2(n); // 2048
  const trueCps = 137.4 / nfft;
  const x = synth(n, trueCps, { dc: 100, ramp: 40, noise: 0.02 });

  const peak = analyseProfile(x);
  assert.ok(peak, 'no peak found');

  const errBins = Math.abs(peak.bin - 137.4);
  assert.ok(errBins < 0.1, `peak off by ${errBins.toFixed(3)} bins (got ${peak.bin.toFixed(3)})`);
  assert.ok(peak.prominence > 5, `prominence too low: ${peak.prominence.toFixed(1)}`);
});

test('interpolation beats raw bin picking on an off-bin tone', () => {
  const n = 1024;
  const nfft = nextPow2(n);
  const trueBin = 80.5; // exactly between two bins, the worst case
  const x = synth(n, trueBin / nfft);
  const { mag, nfft: nf } = spectrum(x);
  const peak = dominantPeak(mag, nf, 4);

  let rawBest = 4;
  for (let k = 4; k < mag.length - 1; k++) if (mag[k] > mag[rawBest]) rawBest = k;

  const rawErr = Math.abs(rawBest - trueBin);
  const interpErr = Math.abs(peak.bin - trueBin);
  assert.ok(interpErr <= rawErr, `interpolation made it worse: ${interpErr} vs ${rawErr}`);
});

test('a flat profile with no tone does not manufacture a confident peak', () => {
  const n = 1080;
  const x = synth(n, 0, { amp: 0, dc: 128, noise: 0.5, seed: 7 });
  const peak = analyseProfile(x);
  // A peak object may still come back; what must not happen is high confidence.
  if (peak) assert.ok(peak.prominence < 5, `manufactured a peak at prominence ${peak.prominence.toFixed(1)}`);
});

test('calibration round-trips: known source in, correct unknown out', () => {
  // Pretend the sensor reads a row every 20 microseconds (50 kHz row rate).
  const trueLineTime = 20e-6;

  // A 120 Hz mains lamp therefore bands at:
  const calCyclesPerRow = 120 * trueLineTime; // 0.0024 cycles/row

  const solved = lineTimeFromKnownSource(calCyclesPerRow, 120);
  assert.ok(Math.abs(solved - trueLineTime) < 1e-12, `line time off: ${solved}`);

  // Now an unknown 400 Hz source through the same sensor.
  const unknownCyclesPerRow = 400 * trueLineTime;
  const hz = hzFromCyclesPerRow(unknownCyclesPerRow, solved);
  assert.ok(Math.abs(hz - 400) < 1e-9, `recovered ${hz} Hz, expected 400`);
});

test('calibration refuses impossible inputs instead of returning a number', () => {
  assert.throws(() => lineTimeFromKnownSource(0.002, 0), /positive/);
  assert.throws(() => lineTimeFromKnownSource(0, 120), /positive/);
  assert.throws(() => hzFromCyclesPerRow(0.002, 0), /positive/);
});

test('end to end on a synthetic 400 Hz scene, recovered within 1 percent', () => {
  const trueLineTime = 20e-6;
  const trueHz = 400;
  const rows = 1080;

  // Build the row profile a camera would actually produce: scene content,
  // vignetting, sensor noise, and the illumination modulation on top.
  const cps = trueHz * trueLineTime;
  const profile = synth(rows, cps, { amp: 3.5, dc: 110, ramp: -35, noise: 1.2, seed: 42 });

  const peak = analyseProfile(profile);
  assert.ok(peak, 'no peak');
  const hz = hzFromCyclesPerRow(peak.cyclesPerSample, trueLineTime);

  const pctErr = (Math.abs(hz - trueHz) / trueHz) * 100;
  assert.ok(pctErr < 1, `recovered ${hz.toFixed(1)} Hz, ${pctErr.toFixed(2)} percent error`);
});

test('aliasNote reports the fold point and the ambiguous twin', () => {
  const lineTime = 20e-6;        // 50 kHz row rate
  const note = aliasNote(400, lineTime);
  assert.ok(Math.abs(note.rowRateHz - 50000) < 1e-6);
  assert.ok(Math.abs(note.nyquistHz - 25000) < 1e-6);
  assert.equal(note.unambiguous, true);
  assert.ok(note.aliases.some((f) => Math.abs(f - 49600) < 1e-6), 'missing the lower alias');

  const folded = aliasNote(30000, lineTime);
  assert.equal(folded.unambiguous, false, 'a frequency above Nyquist must be flagged');
});

// ---------------------------------------------------------------------------
// REGRESSION TESTS from a real false green.
//
// On the first live run against a laptop webcam pointed at an ordinary room,
// the instrument reported "peak / median 394.7", which reads as overwhelming
// confidence, while the row profile was a flat line with one edge in it. Two
// separate defects produced that: a confidence metric that compared low
// frequency scene content against the high-bin noise floor, and a detrend
// window whose corner sat ABOVE the frequency band we were trying to keep.
// ---------------------------------------------------------------------------

/** A scene with no periodic content: smooth gradient plus one sharp edge. */
function sceneProfile(n, { edgeAt = 0.5, noise = 0.6, seed = 3 } = {}) {
  let s = seed >>> 0;
  const rnd = () => ((s = (s * 1664525 + 1013904223) >>> 0) / 4294967296) - 0.5;
  const x = new Float64Array(n);
  const edge = Math.floor(n * edgeAt);
  for (let i = 0; i < n; i++) {
    // A bright upper half falling to a dark lower half, plus vignetting.
    const gradient = 150 - 60 * (i / n) - 25 * Math.sin((Math.PI * i) / n);
    x[i] = gradient + (i > edge ? -40 : 0) + noise * rnd();
  }
  return x;
}

test('REGRESSION: a scene with no periodic content must not score as confident', () => {
  const prof = sceneProfile(1080);
  const peak = analyseProfile(prof);
  if (peak) {
    assert.ok(
      peak.prominence < 8,
      `a blank scene produced prominence ${peak.prominence.toFixed(1)}; this is the 394.7 bug`,
    );
  }
});

test('REGRESSION: the old global-median metric is the thing that was broken', () => {
  // Demonstrates the defect rather than asserting around it. The same scene
  // scores wildly on the global-median ratio and sanely on local prominence.
  const prof = sceneProfile(1080);
  const d = hannCopy(detrend(prof, 401));
  const { mag, nfft } = spectrum(d);

  const globalMetric = dominantPeak(mag, nfft, 3);
  const localMetric = prominentPeak(mag, nfft, { loBin: 3 });

  assert.ok(globalMetric && localMetric);
  assert.ok(
    localMetric.prominence < globalMetric.snr,
    `local prominence (${localMetric.prominence.toFixed(1)}) should be far below the global ratio (${globalMetric.snr.toFixed(1)})`,
  );
});

function hannCopy(x) {
  const y = Float64Array.from(x);
  const n = y.length;
  for (let i = 0; i < n; i++) y[i] *= 0.5 * (1 - Math.cos((2 * Math.PI * i) / (n - 1)));
  return y;
}

// 120 Hz through a sensor reading a row every 27 microseconds, on 1080 rows.
const LINE_TIME = 27e-6;
const MAINS_CPS = 120 * LINE_TIME;   // about 3.24e-3 cycles per row

test('a clean 120 Hz line over a smooth scene is found by the single-strip path', () => {
  const n = 1080;
  // A lit wall: smooth gradient and vignetting, no hard edge.
  const prof = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const gradient = 150 - 40 * (i / n) - 20 * Math.sin((Math.PI * i) / n);
    prof[i] = gradient + 6 * Math.sin(2 * Math.PI * MAINS_CPS * i + 0.7);
  }
  const peak = analyseProfile(prof);
  assert.ok(peak, 'flicker line not found at all');
  const hz = hzFromCyclesPerRow(peak.cyclesPerSample, LINE_TIME);
  assert.ok(Math.abs(hz - 120) < 8, `recovered ${hz.toFixed(1)} Hz, expected 120`);
});

test('DOCUMENTED LIMIT: a hard scene edge defeats the single-strip path at 120 Hz', () => {
  // This is the defect observed on the first real camera run, preserved as a
  // test so nobody "fixes" it by loosening a tolerance. At 1080 rows a 120 Hz
  // lamp completes only ~3.5 cycles, and a step edge carries broadband energy
  // that can out-score it. The instrument must not pretend otherwise, and the
  // page states this limit where a user can read it.
  const n = 1080;
  const scene = sceneProfile(n, { noise: 0.6 });   // contains a step edge
  const prof = new Float64Array(n);
  for (let i = 0; i < n; i++) prof[i] = scene[i] + 6 * Math.sin(2 * Math.PI * MAINS_CPS * i + 0.7);

  const peak = analyseProfile(prof);
  const hz = peak ? hzFromCyclesPerRow(peak.cyclesPerSample, LINE_TIME) : null;
  assert.ok(
    hz === null || Math.abs(hz - 120) > 8,
    `the single-strip path unexpectedly succeeded (${hz?.toFixed(1)} Hz). If this now works, the limit note on the page must be updated.`,
  );

  assert.ok(cyclesInRecord(MAINS_CPS, n) < 4, 'the premise of this limit is fewer than 4 cycles in the record');
});

test('cross-strip coherence recovers the SAME 120 Hz line the single strip lost', () => {
  // Same scene, same lamp, but now we keep the columns separate. The lamp lights
  // every column identically; the edge does not span them.
  const n = 1080;
  const K = 8;
  const profiles = [];
  for (let s = 0; s < K; s++) {
    const scene = sceneProfile(n, { edgeAt: 0.35 + 0.05 * s, noise: 0.6, seed: 11 + s });
    const p = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      // The flicker term is identical in every strip: same amplitude, same phase.
      p[i] = scene[i] + 6 * Math.sin(2 * Math.PI * MAINS_CPS * i + 0.7);
    }
    profiles.push(p);
  }

  const peak = coherentPeak(profiles);
  assert.ok(peak, 'coherent path found nothing');
  const hz = hzFromCyclesPerRow(peak.cyclesPerSample, LINE_TIME);
  assert.ok(Math.abs(hz - 120) < 8, `recovered ${hz.toFixed(1)} Hz, expected 120`);
  assert.ok(peak.coherence > 0.75, `coherence only ${peak.coherence.toFixed(2)}`);
  assert.equal(peak.strips, K);
});

test('coherence is high for a global signal and low for independent phases', () => {
  const n = 512;
  const cps = 0.02;
  const bin = Math.round(cps * nextPow2(n));

  const inPhase = [];
  const randomPhase = [];
  for (let s = 0; s < 8; s++) {
    const a = new Float64Array(n);
    const b = new Float64Array(n);
    const phi = s * 1.31;
    for (let i = 0; i < n; i++) {
      a[i] = Math.sin(2 * Math.PI * cps * i);        // same phase everywhere
      b[i] = Math.sin(2 * Math.PI * cps * i + phi);  // scattered phases
    }
    inPhase.push(spectrumComplex(a));
    randomPhase.push(spectrumComplex(b));
  }

  const cGlobal = coherence(inPhase, bin);
  const cLocal = coherence(randomPhase, bin);
  assert.ok(cGlobal > 0.98, `a global signal should be near 1, got ${cGlobal.toFixed(3)}`);
  assert.ok(cLocal < 0.6, `scattered phases should be low, got ${cLocal.toFixed(3)}`);
});

test('coherentPeak refuses a scene with no global periodic content', () => {
  const n = 1080;
  const profiles = [];
  for (let s = 0; s < 8; s++) profiles.push(sceneProfile(n, { edgeAt: 0.3 + 0.06 * s, seed: 40 + s }));
  const peak = coherentPeak(profiles, { minCoherence: 0.75 });
  if (peak) {
    assert.ok(peak.coherence < 0.95, `pure scenery reported coherence ${peak.coherence.toFixed(2)}`);
  }
});

test('cyclesInRecord reports how much the instrument actually saw', () => {
  assert.ok(Math.abs(cyclesInRecord(MAINS_CPS, 1080) - 3.5) < 0.2);
  assert.ok(cyclesInRecord(MAINS_CPS, 1080 * 8) > 25, 'stacking frames buys cycles');
});

test('REGRESSION: the old 65-row detrend window destroys the 120 Hz band', () => {
  // This is the second defect, stated as an executable fact. A 65-row window
  // has its corner near 1/65 cycles per row, which is well above the roughly
  // 3.2e-3 cycles per row where a mains lamp actually lives.
  const lineTime = 27e-6;
  const cps = 120 * lineTime;
  const n = 1080;
  const prof = new Float64Array(n);
  for (let i = 0; i < n; i++) prof[i] = 120 + 6 * Math.sin(2 * Math.PI * cps * i);

  const keptShort = spectrum(hannCopy(detrend(prof, 65)));
  const keptLong = spectrum(hannCopy(detrend(prof, 401)));
  const bin = Math.round(cps * keptLong.nfft);

  assert.ok(
    keptLong.mag[bin] > keptShort.mag[bin] * 3,
    `the long window must preserve far more of the signal (${keptLong.mag[bin].toFixed(1)} vs ${keptShort.mag[bin].toFixed(1)})`,
  );
});

test('detrendWindowFor derives a window from the lowest frequency to keep', () => {
  const cps = 120 * 27e-6;
  const win = detrendWindowFor(cps, 4);
  assert.ok(win > 1000, `window ${win} is too short to pass a 120 Hz lamp`);
  assert.equal(win % 2, 1, 'window should be odd so the average is centred');
  assert.throws(() => detrendWindowFor(0), /positive/);
});

test('localBaseline tracks a hump but does not swallow an isolated line', () => {
  const n = 512;
  const mag = new Float64Array(n);
  for (let k = 0; k < n; k++) mag[k] = 10 * Math.exp(-k / 60) + 0.5; // scene hump
  mag[200] = 40;                                                      // isolated line
  const base = localBaseline(mag, 24);
  assert.ok(base[200] < 5, `baseline absorbed the line: ${base[200]}`);
  assert.ok(base[10] > 3, `baseline failed to follow the hump: ${base[10]}`);
});

test('PeakTracker withholds a verdict until it has repeated', () => {
  const t = new PeakTracker(12);
  assert.equal(t.verdict(), null, 'a verdict with no samples is not a verdict');
  for (let i = 0; i < 3; i++) t.push(0.00324);
  assert.equal(t.verdict(), null, 'three frames is still not enough');
  for (let i = 0; i < 6; i++) t.push(0.00324);
  const v = t.verdict();
  assert.ok(v && v.stable, 'a steady line should be called stable');
  assert.ok(Math.abs(v.median - 0.00324) < 1e-9);
});

test('PeakTracker refuses to call a jittering estimate stable', () => {
  const t = new PeakTracker(12);
  const jumpy = [0.0030, 0.0071, 0.0042, 0.0095, 0.0033, 0.0088, 0.0051, 0.0062];
  jumpy.forEach((v) => t.push(v));
  const v = t.verdict();
  assert.ok(v, 'should return a verdict object');
  assert.equal(v.stable, false, `spread ${v.spreadPct.toFixed(1)}% was wrongly called stable`);
});
