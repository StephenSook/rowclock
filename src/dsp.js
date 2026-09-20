// dsp.js - the numerical core. No dependencies, no build step, no network.
//
// Everything here operates on plain Float64Array. It is written to be callable
// identically from the browser page and from a Node harness, so the same code
// produces the numbers in the live demo and the numbers in the tests.

/**
 * In-place iterative radix-2 Cooley-Tukey FFT.
 * re and im must be the same length and that length must be a power of two.
 * Decimation in time, with the bit-reversal permutation done up front.
 */
export function fft(re, im) {
  const n = re.length;
  if (n !== im.length) throw new Error(`fft: re/im length mismatch (${n} vs ${im.length})`);
  if (n < 2 || (n & (n - 1)) !== 0) throw new Error(`fft: length must be a power of two, got ${n}`);

  // Bit-reversal permutation.
  for (let i = 1, j = 0; i < n; i++) {
    let bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      let t = re[i]; re[i] = re[j]; re[j] = t;
      t = im[i]; im[i] = im[j]; im[j] = t;
    }
  }

  // Butterflies.
  for (let len = 2; len <= n; len <<= 1) {
    const ang = (-2 * Math.PI) / len;
    const wRe = Math.cos(ang);
    const wIm = Math.sin(ang);
    for (let i = 0; i < n; i += len) {
      let curRe = 1, curIm = 0;
      for (let k = 0; k < len / 2; k++) {
        const aRe = re[i + k],            aIm = im[i + k];
        const bRe = re[i + k + len / 2],  bIm = im[i + k + len / 2];
        const tRe = bRe * curRe - bIm * curIm;
        const tIm = bRe * curIm + bIm * curRe;
        re[i + k] = aRe + tRe;  im[i + k] = aIm + tIm;
        re[i + k + len / 2] = aRe - tRe;
        im[i + k + len / 2] = aIm - tIm;
        const nextRe = curRe * wRe - curIm * wIm;
        curIm = curRe * wIm + curIm * wRe;
        curRe = nextRe;
      }
    }
  }
}

/**
 * Inverse FFT, in place, by conjugation.
 *
 * ifft(X) = conj(fft(conj(X))) / n. Reusing the forward transform means there
 * is one butterfly implementation to get right and to test, rather than two.
 */
export function ifft(re, im) {
  const n = re.length;
  for (let i = 0; i < n; i++) im[i] = -im[i];
  fft(re, im);
  for (let i = 0; i < n; i++) {
    re[i] /= n;
    im[i] = -im[i] / n;
  }
}

/** Next power of two >= n. */
export function nextPow2(n) {
  let p = 1;
  while (p < n) p <<= 1;
  return p;
}

/**
 * Remove DC and low-order scene structure.
 *
 * This matters more than the FFT does. A row-mean profile of a real photograph
 * is dominated by the SCENE (a bright ceiling, a dark desk, lens vignetting),
 * and that content sits at very low spatial frequency. Subtracting a centered
 * moving average of width `win` is a crude high-pass that leaves the flicker
 * band untouched while removing the picture.
 *
 * Returns a new Float64Array; does not modify the input.
 */
export function detrend(x, win = 65) {
  const n = x.length;
  const out = new Float64Array(n);
  if (win < 3) { out.set(x); return out; }
  const half = win >> 1;

  // Prefix sums so the moving average is O(n) rather than O(n*win).
  const pre = new Float64Array(n + 1);
  for (let i = 0; i < n; i++) pre[i + 1] = pre[i] + x[i];

  for (let i = 0; i < n; i++) {
    const a = Math.max(0, i - half);
    const b = Math.min(n, i + half + 1);
    const mean = (pre[b] - pre[a]) / (b - a);
    out[i] = x[i] - mean;
  }
  return out;
}

/** Hann window, applied in place. Reduces spectral leakage from the finite record. */
export function hann(x) {
  const n = x.length;
  for (let i = 0; i < n; i++) {
    x[i] *= 0.5 * (1 - Math.cos((2 * Math.PI * i) / (n - 1)));
  }
  return x;
}

/**
 * Magnitude spectrum of a real signal, zero-padded to the next power of two.
 * Returns { mag, nfft } where mag has nfft/2 usable bins.
 * Bin k corresponds to k / nfft cycles per sample.
 */
export function spectrum(x) {
  const nfft = nextPow2(x.length);
  const re = new Float64Array(nfft);
  const im = new Float64Array(nfft);
  re.set(x);
  fft(re, im);
  const half = nfft >> 1;
  const mag = new Float64Array(half);
  for (let k = 0; k < half; k++) {
    mag[k] = Math.hypot(re[k], im[k]);
  }
  return { mag, nfft };
}

/**
 * Running-median baseline of a spectrum.
 *
 * This exists because comparing a peak to the GLOBAL median is not a test of
 * anything. A real scene's row profile is dominated by low spatial frequency
 * (a bright ceiling, a dark desk, a horizon line), so the low bins tower over
 * the high-bin noise floor and the global ratio comes out enormous no matter
 * what the camera is pointed at. Measured on a blank wall it read 394.7, which
 * looks like overwhelming confidence and means nothing.
 *
 * A LOCAL baseline asks the only question that matters: does this bin stand
 * above its own neighbourhood. A periodic source produces an isolated line
 * that does. Smooth scene structure does not, because its neighbours are just
 * as large.
 */
export function localBaseline(mag, halfWidth = 24, guard = 2) {
  const n = mag.length;
  const base = new Float64Array(n);
  const buf = [];
  for (let k = 0; k < n; k++) {
    // The window narrows near DC. With a fixed wide window a genuine line at
    // bin 7 has the whole steep scene hump inside its own baseline, which
    // inflates the baseline and buries the line. Measured consequence: a real
    // 120 Hz line at bin ~7 lost to a 4-row pipeline artifact at bin 512.
    const hw = Math.max(6, Math.min(halfWidth, Math.round(k * 0.6)));
    const a = Math.max(0, k - hw);
    const b = Math.min(n, k + hw + 1);
    buf.length = 0;
    for (let j = a; j < b; j++) {
      // A guard band, so a line does not raise the baseline it is measured
      // against. Without it a strong narrow peak partly hides itself.
      if (Math.abs(j - k) <= guard) continue;
      buf.push(mag[j]);
    }
    if (!buf.length) { base[k] = mag[k]; continue; }
    buf.sort((p, q) => p - q);
    base[k] = buf[buf.length >> 1];
  }
  return base;
}

/**
 * Locate the dominant peak and refine it to sub-bin precision.
 *
 * Quadratic interpolation over the log-magnitudes of the peak bin and its two
 * neighbours. Without this the frequency resolution is one bin, which at
 * nfft=2048 over 1080 rows is a coarse answer; with it the estimate is good to
 * a small fraction of a bin, which is what makes a three-significant-figure
 * readout honest rather than decorative.
 *
 * `loBin` skips the near-DC bins that survive detrending.
 * Returns { bin, cyclesPerSample, magnitude, snr } or null if nothing stands out.
 */
export function dominantPeak(mag, nfft, loBin = 4) {
  const n = mag.length;
  if (loBin >= n - 1) return null;

  let best = loBin;
  for (let k = loBin; k < n - 1; k++) {
    if (mag[k] > mag[best]) best = k;
  }
  if (best <= 0 || best >= n - 1) return null;
  if (!(mag[best] > 0)) return null;

  // Quadratic (parabolic) interpolation in the log domain.
  const eps = 1e-12;
  const a = Math.log(mag[best - 1] + eps);
  const b = Math.log(mag[best] + eps);
  const c = Math.log(mag[best + 1] + eps);
  const denom = a - 2 * b + c;
  const delta = denom === 0 ? 0 : (0.5 * (a - c)) / denom;
  const refined = best + Math.max(-0.5, Math.min(0.5, delta));

  // Median of the whole band as a noise floor, so the ratio is robust to the
  // peak itself and to a handful of other strong bins.
  const sorted = Float64Array.from(mag.subarray(loBin)).sort();
  const median = sorted[sorted.length >> 1] || eps;

  return {
    bin: refined,
    cyclesPerSample: refined / nfft,
    magnitude: mag[best],
    snr: mag[best] / median,
  };
}

/**
 * Find the most PROMINENT line in the spectrum, not merely the tallest bin.
 *
 * Prominence is magnitude divided by the local baseline. An isolated periodic
 * line scores high; the shoulder of a broad scene-content hump scores near 1
 * however tall it is in absolute terms. This is the honest confidence number
 * and it is what gates whether a reading is shown at all.
 */
export function prominentPeak(mag, nfft, { loBin = 3, baselineHalfWidth = 24 } = {}) {
  const n = mag.length;
  if (loBin >= n - 1) return null;
  const base = localBaseline(mag, baselineHalfWidth);

  let best = -1;
  let bestRatio = 0;
  for (let k = loBin; k < n - 1; k++) {
    // Only consider actual local maxima, so we never report the flank of a hump.
    if (!(mag[k] >= mag[k - 1] && mag[k] >= mag[k + 1])) continue;
    const ratio = mag[k] / (base[k] + 1e-12);
    if (ratio > bestRatio) { bestRatio = ratio; best = k; }
  }
  if (best < 1 || best >= n - 1) return null;

  const eps = 1e-12;
  const a = Math.log(mag[best - 1] + eps);
  const b = Math.log(mag[best] + eps);
  const c = Math.log(mag[best + 1] + eps);
  const denom = a - 2 * b + c;
  const delta = denom === 0 ? 0 : (0.5 * (a - c)) / denom;
  const refined = best + Math.max(-0.5, Math.min(0.5, delta));

  return {
    bin: refined,
    cyclesPerSample: refined / nfft,
    magnitude: mag[best],
    prominence: bestRatio,
  };
}

/** Complex spectrum of a real signal, zero-padded. Returns { re, im, nfft }. */
export function spectrumComplex(x) {
  const nfft = nextPow2(x.length);
  const re = new Float64Array(nfft);
  const im = new Float64Array(nfft);
  re.set(x);
  fft(re, im);
  return { re, im, nfft };
}

/**
 * CROSS-STRIP PHASE COHERENCE, which is the discriminator that makes a single
 * frame usable at low cycle counts.
 *
 * The problem this solves: at 1080 rows and a typical row time, a 120 Hz mains
 * lamp completes only about three and a half cycles across the whole frame.
 * That is barely enough to call periodic, and a single hard edge in the scene
 * (a wall meeting a ceiling) carries broadband energy that can out-score it.
 * Measured on the first real build, exactly that happened: a 120 Hz source was
 * reported as 71.7 Hz because the picker preferred the edge.
 *
 * The way out is physical rather than numerical. Illumination flicker is a
 * GLOBAL property of the frame: every column is lit by the same lamp, so the
 * banding has the same phase everywhere across the width. Scene structure is
 * LOCAL: an edge occupies some columns and not others, and its phase varies
 * across the width or is absent entirely.
 *
 * So: split the frame into K vertical strips, transform each strip's own row
 * profile, and at every candidate bin measure how well the K complex values
 * agree in phase. This is the phase-locking value, |sum(Z)| / sum(|Z|), which
 * is 1 when all strips agree exactly and falls toward 1/sqrt(K) for
 * independent phases.
 *
 * A candidate now has to be BOTH prominent above its local baseline AND
 * coherent across the frame's width. An edge fails the second test even when
 * it wins the first.
 */
export function coherence(strips, bin) {
  const k = Math.round(bin);
  let sumRe = 0, sumIm = 0, sumMag = 0;
  for (const s of strips) {
    if (k >= s.re.length) return 0;
    const re = s.re[k], im = s.im[k];
    sumRe += re; sumIm += im;
    sumMag += Math.hypot(re, im);
  }
  if (sumMag <= 0) return 0;
  return Math.hypot(sumRe, sumIm) / sumMag;
}

/**
 * Find the line that is both prominent and coherent across the frame width.
 *
 * `profiles` is an array of K per-strip row profiles, already the same length.
 * Returns the winning candidate with both scores attached, or null.
 */
export function coherentCandidates(profiles, {
  loBin = 3,
  hiCyclesPerRow = 0.2,
  baselineHalfWidth = 24,
  detrendWin = 401,
  minCoherence = 0.75,
  minProminence = 2,
  limit = 6,
} = {}) {
  if (!profiles?.length) return [];

  const strips = profiles.map((p) => spectrumComplex(hann(detrend(p, detrendWin))));
  const nfft = strips[0].nfft;
  const half = nfft >> 1;

  // The full-width profile is the mean of the strips, and it is what we score
  // prominence on, because averaging is what buys the signal-to-noise.
  const mean = new Float64Array(half);
  for (let k = 0; k < half; k++) {
    let acc = 0;
    for (const s of strips) acc += Math.hypot(s.re[k], s.im[k]);
    mean[k] = acc / strips.length;
  }

  const base = localBaseline(mean, baselineHalfWidth);

  // An upper bound in cycles PER ROW, not in hertz, because it is a statement
  // about the row grid rather than about the world. Above roughly 0.2 cycles
  // per row a "line" is a pattern repeating every few rows, which is the
  // signature of the imaging pipeline itself (demosaic, chroma upsampling,
  // rescaling) and not of a lamp. Measured case: a 4-row period at bin 512
  // with coherence 0.947 and frame-to-frame spread 0.04 percent, which is
  // exactly how a FIXED spatial pattern behaves, since it is perfectly stable
  // and perfectly coherent by construction.
  const hiBin = Math.min(half - 2, Math.floor(hiCyclesPerRow * nfft));

  const out = [];
  for (let k = loBin; k < hiBin; k++) {
    if (!(mean[k] >= mean[k - 1] && mean[k] >= mean[k + 1])) continue;
    const prom = mean[k] / (base[k] + 1e-12);
    if (prom < minProminence) continue;
    const coh = coherence(strips, k);
    if (coh < minCoherence) continue;

    const eps = 1e-12;
    const a = Math.log(mean[k - 1] + eps);
    const b = Math.log(mean[k] + eps);
    const c = Math.log(mean[k + 1] + eps);
    const denom = a - 2 * b + c;
    const delta = denom === 0 ? 0 : (0.5 * (a - c)) / denom;
    const refined = k + Math.max(-0.5, Math.min(0.5, delta));

    out.push({
      bin: refined,
      cyclesPerSample: refined / nfft,
      prominence: prom,
      coherence: coh,
      score: prom * coh,
      strips: strips.length,
      nfft,
    });
  }

  out.sort((p, q) => q.score - p.score);
  return out.slice(0, limit);
}

/** The single best candidate, or null. */
export function coherentPeak(profiles, opts = {}) {
  const c = coherentCandidates(profiles, { ...opts, limit: 1 });
  return c.length ? c[0] : null;
}

/**
 * Decide whether a candidate is a temporal signal or a fixed spatial pattern.
 *
 * This is the distinguishing-prediction test, and it is the only honest way to
 * tell the two apart from inside a browser.
 *
 * A TEMPORAL signal at frequency f produces cyclesPerRow = f * lineTime. The
 * row time depends on the readout mode, so changing the capture resolution
 * changes lineTime, and therefore changes cyclesPerRow, while f stays put.
 *
 * A SPATIAL artifact is locked to the row grid. It repeats every N rows
 * whatever the resolution, so its cyclesPerRow does not move.
 *
 * Measure the same scene at two resolutions and the two hypotheses predict
 * different things. That is a real experiment, not an assertion, and a reading
 * that has not passed it should not be presented as a measurement.
 */
export function classifyByResolutionSwap(cyclesPerRowA, cyclesPerRowB, { tolPct = 5 } = {}) {
  if (!(cyclesPerRowA > 0) || !(cyclesPerRowB > 0)) {
    return { verdict: 'inconclusive', reason: 'one of the two captures produced no stable line' };
  }
  const changePct = (Math.abs(cyclesPerRowB - cyclesPerRowA) / cyclesPerRowA) * 100;
  if (changePct <= tolPct) {
    return {
      verdict: 'spatial',
      changePct,
      reason: `cycles per row barely moved (${changePct.toFixed(1)}%) when the row grid changed. This is locked to the sensor grid, so it is a pipeline artifact, not a light.`,
    };
  }
  return {
    verdict: 'temporal',
    changePct,
    reason: `cycles per row moved ${changePct.toFixed(1)}% with the row time, which is what a real temporal signal does and a fixed pattern cannot.`,
  };
}

/**
 * Reconstruct the sinusoid the detected bin actually represents, and report
 * how much of the signal it explains.
 *
 * This exists because every confidence number we display can be fooled by the
 * same thing. A genuine 120 Hz banding pattern and a single drifting
 * horizontal edge both span the full width, so both score high coherence;
 * both are stable frame to frame; and an isolated bump roughly one sixth of
 * the frame tall puts its spectral peak near bin 6 for reasons that have
 * nothing to do with a lamp. Observed on a real ceiling: bin 6.37, coherence
 * 0.970, spread 1.31 percent, and a row profile showing ONE bump rather than
 * the three ripples that frequency demands.
 *
 * A number cannot settle that. Drawing the fitted wave over the data can: if
 * the line is real the wave tracks the whole trace, and if it is a lump the
 * wave oscillates through flat regions where nothing is happening.
 *
 * Returns the reconstruction plus fitFraction, the share of the signal's
 * energy the single sinusoid accounts for. A periodic line explains most of
 * it; a localized bump explains very little, because its energy is spread
 * across many bins.
 */
export function fitSinusoid(signal, bin, nfft) {
  const n = signal.length;
  const w = (2 * Math.PI * bin) / nfft;

  // Least squares at the EXACT frequency, not the value of one rounded FFT
  // bin. The detected bin is deliberately fractional, and a fractional tone
  // leaks across neighbours, so a single bin's magnitude understates the
  // amplitude badly: a pure line fitted that way explained only 65 percent of
  // its own energy. Projecting onto cos and sin at the true frequency is
  // exact, costs one pass, and does not care about zero padding.
  let cc = 0, ss = 0, cs = 0, xc = 0, xs = 0;
  for (let i = 0; i < n; i++) {
    const c = Math.cos(w * i);
    const s = Math.sin(w * i);
    cc += c * c; ss += s * s; cs += c * s;
    xc += signal[i] * c; xs += signal[i] * s;
  }
  const det = cc * ss - cs * cs;
  if (Math.abs(det) < 1e-12) {
    return { fitted: new Float64Array(n), amplitude: 0, phase: 0, fitFraction: 0 };
  }
  const A = (xc * ss - xs * cs) / det;
  const B = (xs * cc - xc * cs) / det;

  const amp = Math.hypot(A, B);
  const phase = Math.atan2(-B, A);

  const fitted = new Float64Array(n);
  for (let i = 0; i < n; i++) fitted[i] = A * Math.cos(w * i) + B * Math.sin(w * i);

  let sigE = 0, resE = 0;
  for (let i = 0; i < n; i++) {
    sigE += signal[i] * signal[i];
    const d = signal[i] - fitted[i];
    resE += d * d;
  }
  const fitFraction = sigE > 0 ? Math.max(0, 1 - resE / sigE) : 0;

  return { fitted, amplitude: amp, phase, fitFraction };
}

/**
 * Running mean of the row profile, which IS the static scene.
 *
 * This is the move that makes the instrument work, and it took two real camera
 * runs to find. The diagnostic that produced it: a genuine 120 Hz line scored
 * prominence 1.77 and coherence 0.703, both below bar, because with only about
 * three and a half cycles across 1080 rows the line is intrinsically four bins
 * wide and sits on top of strong scene energy. Tuning thresholds cannot
 * separate them; they overlap in frequency.
 *
 * What separates them is TIME, not frequency:
 *
 *   - The scene is static. A wall, an edge, vignetting and the sensor's own
 *     fixed pattern produce the same row profile every frame, with the same
 *     phase.
 *   - Flicker is not. The frame period is not locked to the mains, so the
 *     banding sits in a different place each frame and its phase walks.
 *
 * So the running mean converges on the scene and the flicker averages itself
 * away. Subtract it, and what remains is the part of the image that MOVES.
 * That removes the gradient, the edges, the vignetting, and, for free, the
 * four-row pipeline artifact that nearly became our headline number, because
 * a fixed pattern is static by definition.
 *
 * Camera motion breaks the static assumption, which is why motionEnergy() is
 * published alongside the reading rather than hidden.
 */
export class TemporalBackground {
  constructor(alpha = 0.06) {
    this.alpha = alpha;   // exponential forgetting factor
    this.mean = null;
    this.frames = 0;
  }

  reset() { this.mean = null; this.frames = 0; }

  /** Update with a new profile and return the residual (profile minus scene). */
  update(profile) {
    const n = profile.length;
    if (!this.mean || this.mean.length !== n) {
      this.mean = Float64Array.from(profile);
      this.frames = 1;
      return new Float64Array(n); // nothing to say from a single frame
    }
    const a = this.alpha;
    const resid = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      resid[i] = profile[i] - this.mean[i];
      this.mean[i] = (1 - a) * this.mean[i] + a * profile[i];
    }
    this.frames++;
    return resid;
  }

  /** True once the running mean has seen enough frames to be a scene estimate. */
  get ready() { return this.frames >= 8; }
}

/**
 * How much of the residual is bulk movement rather than modulation.
 *
 * If the camera is being waved around, every row changes and the residual is
 * large and broadband. A reading taken then is not measuring a lamp, it is
 * measuring a shaky hand, and the page says so instead of printing a number.
 */
export function motionEnergy(residual, profile) {
  let rs = 0, ps = 0;
  for (let i = 0; i < residual.length; i++) { rs += residual[i] * residual[i]; ps += profile[i] * profile[i]; }
  return ps > 0 ? Math.sqrt(rs / ps) : 0;
}

/**
 * How many cycles of a given frequency fit in the record.
 *
 * Published beside every reading, because it is the honest statement of how
 * much the instrument actually saw. Below about four cycles a periodogram
 * estimate is weakly determined no matter how clean the arithmetic is, and a
 * reading that does not disclose this is overclaiming.
 */
export function cyclesInRecord(cyclesPerSample, nSamples) {
  return cyclesPerSample * nSamples;
}

/**
 * Choose a detrend window from the LOWEST frequency we intend to keep.
 *
 * Getting this backwards is how the first build of this instrument deleted its
 * own signal. A moving-average high-pass has its corner near one over the
 * window length, so a 65-row window rejects everything below about 1/65
 * cycles per row. A 120 Hz lamp at roughly 27 microseconds per row bands at
 * about 3.2e-3 cycles per row, which is a factor of five BELOW that corner.
 * The filter was sitting on top of the thing it was supposed to pass.
 *
 * The window must span several periods of the slowest signal of interest.
 */
export function detrendWindowFor(lowestCyclesPerRow, periods = 4) {
  if (!(lowestCyclesPerRow > 0)) throw new Error('detrendWindowFor: need a positive cycles/row');
  const period = 1 / lowestCyclesPerRow;
  const win = Math.round(period * periods);
  return Math.max(3, win | 1); // odd, so the moving average is centred
}

/**
 * The whole pipeline for one profile: detrend, window, transform, pick the peak.
 *
 * `detrendWin` defaults to 401 rows, which passes everything above roughly
 * 2.5e-3 cycles per row and therefore keeps mains flicker on a 1080-row frame,
 * while still removing gross scene gradient and vignetting.
 */
export function analyseProfile(profile, { detrendWin = 401, loBin = 3, baselineHalfWidth = 24 } = {}) {
  if (!profile || profile.length < 32) return null;
  const d = hann(detrend(profile, detrendWin));
  const { mag, nfft } = spectrum(d);
  return prominentPeak(mag, nfft, { loBin, baselineHalfWidth });
}

/**
 * Cross-frame agreement, which is the discriminator that actually works.
 *
 * Mains flicker is stationary: the same line appears frame after frame in the
 * same place. Scene content is not, because the camera shakes, the subject
 * moves, and auto-gain breathes. So a reading earns trust by REPEATING, and a
 * single frame never earns it at all.
 *
 * Reports the median estimate and the spread as a percentage of it.
 */
export class PeakTracker {
  constructor(depth = 12) {
    this.depth = depth;
    this.samples = [];
  }

  push(cyclesPerSample) {
    if (!(cyclesPerSample > 0)) return;
    this.samples.push(cyclesPerSample);
    if (this.samples.length > this.depth) this.samples.shift();
  }

  reset() { this.samples = []; }

  /** { median, spreadPct, n, stable } or null while still filling. */
  verdict({ minSamples = 6, maxSpreadPct = 2 } = {}) {
    const n = this.samples.length;
    if (n < minSamples) return null;
    const s = [...this.samples].sort((a, b) => a - b);
    const median = s[n >> 1];
    // Interquartile spread, so one bad frame does not veto a good run.
    const q1 = s[Math.floor(n * 0.25)];
    const q3 = s[Math.floor(n * 0.75)];
    const spreadPct = median > 0 ? ((q3 - q1) / median) * 100 : Infinity;
    return { median, spreadPct, n, stable: spreadPct <= maxSpreadPct };
  }
}

/**
 * Solve the sensor's per-row readout time from a source of KNOWN frequency.
 *
 * This is the calibration that makes the whole instrument possible, and it is
 * why no datasheet is needed. A mains-powered lamp modulates at twice the line
 * frequency (120 Hz on a 60 Hz supply) because instantaneous power goes as
 * v(t)^2 and therefore peaks twice per cycle. Photograph one, measure the
 * banding period in CYCLES PER ROW, and the row period falls out:
 *
 *     cyclesPerRow [cycles/row] = f_known [cycles/s] * lineTime [s/row]
 *  => lineTime = cyclesPerRow / f_known
 *
 * Returns seconds per row.
 */
export function lineTimeFromKnownSource(cyclesPerRow, knownHz) {
  if (!(knownHz > 0)) throw new Error('lineTimeFromKnownSource: knownHz must be positive');
  if (!(cyclesPerRow > 0)) throw new Error('lineTimeFromKnownSource: cyclesPerRow must be positive');
  return cyclesPerRow / knownHz;
}

/** Convert a measured spatial frequency to a temporal one, given the calibration. */
export function hzFromCyclesPerRow(cyclesPerRow, lineTimeSeconds) {
  if (!(lineTimeSeconds > 0)) throw new Error('hzFromCyclesPerRow: lineTime must be positive');
  return cyclesPerRow / lineTimeSeconds;
}

/**
 * The frequency the rolling shutter CANNOT distinguish from f.
 *
 * Row sampling is still sampling, so it aliases. With a row rate of
 * 1/lineTime, anything above half that folds back. Reporting the fold points
 * alongside a reading is the difference between an instrument and a number
 * generator, and it is the first thing an imaging engineer will ask about.
 */
export function aliasNote(hz, lineTimeSeconds) {
  const rowRate = 1 / lineTimeSeconds;
  const nyquist = rowRate / 2;
  return {
    rowRateHz: rowRate,
    nyquistHz: nyquist,
    unambiguous: hz < nyquist,
    // The other candidates that produce an identical banding period.
    aliases: [rowRate - hz, rowRate + hz].filter((f) => f > 0 && f < rowRate * 3),
  };
}
