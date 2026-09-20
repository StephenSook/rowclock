// capture.js - getting an honest row profile out of a browser camera.
//
// The whole instrument depends on one thing being true: that the frames a
// browser hands us still carry the row-by-row timing structure the sensor
// produced. Phone camera stacks fight us here. They apply auto-exposure, they
// apply temporal noise reduction, and several of them apply explicit
// anti-banding (a "flicker correction") whose entire purpose is to remove the
// signal we are trying to measure.
//
// So this module does two jobs: ask the platform, as loudly as the spec allows,
// to stop doing that; and then report truthfully what the platform actually
// agreed to, so a reading is never presented as trustworthy when the capture
// path was compromised.

/**
 * Open a camera stream at the highest row count we can get.
 *
 * More rows is strictly better: rows are our samples, so a 1080-row frame is a
 * 1080-point record and a 480-row frame is a 480-point one.
 */
export async function openCamera({ facingMode = 'environment', height = 1080 } = {}) {
  if (!navigator.mediaDevices?.getUserMedia) {
    throw new Error('This browser has no getUserMedia. Camera capture is unavailable.');
  }

  const constraints = {
    audio: false,
    video: {
      facingMode,
      width: { ideal: 1920 },
      height: { ideal: height },
      frameRate: { ideal: 30 },
    },
  };

  const stream = await navigator.mediaDevices.getUserMedia(constraints);
  const track = stream.getVideoTracks()[0];

  // Best effort: turn off the processing that would erase the signal. These
  // are optional constraints, so a platform that does not support them simply
  // ignores them rather than failing the call. We read back what stuck.
  const requested = {};
  try {
    const caps = track.getCapabilities ? track.getCapabilities() : {};
    const advanced = [];

    if (caps.exposureMode?.includes('manual')) {
      advanced.push({ exposureMode: 'manual' });
      requested.exposureMode = 'manual';
    }
    if (caps.whiteBalanceMode?.includes('manual')) {
      advanced.push({ whiteBalanceMode: 'manual' });
      requested.whiteBalanceMode = 'manual';
    }
    // powerLineFrequency is the anti-banding control. Asking for "none" asks
    // the stack to stop cancelling mains flicker, which is the exact signal we
    // want to keep.
    if (caps.powerLineFrequency) {
      advanced.push({ powerLineFrequency: 0 });
      requested.powerLineFrequency = 0;
    }

    if (advanced.length) await track.applyConstraints({ advanced });
  } catch {
    // Constraint application is advisory. A refusal is information, not a fault.
  }

  return { stream, track, requested };
}

/**
 * What the platform actually gave us, as opposed to what we asked for.
 *
 * This is reported on screen next to every reading. A number produced through
 * a pipeline that is silently cancelling flicker is not a measurement, and the
 * user is entitled to see which of the protections we actually got.
 */
export function captureHealth(track) {
  const settings = track.getSettings ? track.getSettings() : {};
  const caps = track.getCapabilities ? track.getCapabilities() : {};
  const warnings = [];

  if (settings.exposureMode && settings.exposureMode !== 'manual') {
    warnings.push('Auto-exposure is on. The stack may be varying exposure between frames.');
  }
  if (!('exposureMode' in caps)) {
    warnings.push('This platform exposes no exposure control at all.');
  }
  if (settings.powerLineFrequency) {
    warnings.push(
      `Anti-banding is active at ${settings.powerLineFrequency} Hz. The camera is actively cancelling the signal this instrument reads.`,
    );
  }
  if ((settings.height || 0) < 720) {
    warnings.push(`Only ${settings.height} rows per frame. Fewer rows means a coarser spectrum.`);
  }

  return { settings, capabilities: caps, warnings };
}

/**
 * Pull one frame and reduce it to a per-row brightness profile.
 *
 * Green channel, not luminance. On a Bayer sensor green is sampled at twice
 * the rate of red or blue, so it carries the best signal-to-noise per row, and
 * using it avoids paying for a colour conversion on every pixel.
 *
 * `columnFraction` narrows the horizontal span that is averaged. Averaging the
 * full width is right for an illumination measurement, where the modulation is
 * global. It is wrong when the thing of interest occupies part of the frame,
 * because the rest of the scene then dilutes it.
 */
export function rowProfile(ctx, width, height, { columnFraction = 1 } = {}) {
  const span = Math.max(1, Math.round(width * columnFraction));
  const x0 = Math.floor((width - span) / 2);

  const img = ctx.getImageData(x0, 0, span, height);
  const d = img.data;
  const profile = new Float64Array(height);

  for (let y = 0; y < height; y++) {
    let sum = 0;
    const base = y * span * 4;
    for (let x = 0; x < span; x++) {
      sum += d[base + x * 4 + 1]; // green
    }
    profile[y] = sum / span;
  }
  return profile;
}

/**
 * Split the frame into K vertical strips and return a row profile for each.
 *
 * This is what makes the coherence test possible. Illumination flicker lights
 * every column identically, so its banding has the same phase in every strip.
 * Scene structure does not: an edge, a shadow, a dark object occupies some
 * columns and not others. Keeping the strips separate preserves exactly the
 * information that averaging the full width throws away.
 *
 * One pass over the pixels, K accumulators per row, so this costs essentially
 * the same as the single full-width profile it replaces.
 */
export function rowProfileStrips(ctx, width, height, strips = 8) {
  const img = ctx.getImageData(0, 0, width, height);
  const d = img.data;
  const out = [];
  for (let s = 0; s < strips; s++) out.push(new Float64Array(height));

  const stripW = width / strips;

  for (let y = 0; y < height; y++) {
    const rowBase = y * width * 4;
    const sums = new Float64Array(strips);
    const counts = new Int32Array(strips);
    for (let x = 0; x < width; x++) {
      const s = Math.min(strips - 1, (x / stripW) | 0);
      sums[s] += d[rowBase + x * 4 + 1]; // green
      counts[s]++;
    }
    for (let s = 0; s < strips; s++) out[s][y] = counts[s] ? sums[s] / counts[s] : 0;
  }
  return out;
}

/** Mean of a set of strip profiles, which is the full-width profile. */
export function meanProfile(profiles) {
  const n = profiles[0].length;
  const out = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    let acc = 0;
    for (const p of profiles) acc += p[i];
    out[i] = acc / profiles.length;
  }
  return out;
}

/**
 * Saturation check.
 *
 * A clipped row carries no modulation information, because the sensor stopped
 * responding. If a meaningful share of the frame is at the rail, the reading
 * is not merely noisy, it is measuring a flat line, and the honest response is
 * to say so rather than to report the spectrum of the clipping.
 */
export function clippingFraction(ctx, width, height, samples = 4000) {
  const img = ctx.getImageData(0, 0, width, height);
  const d = img.data;
  const total = width * height;
  const stride = Math.max(1, Math.floor(total / samples));
  let clipped = 0;
  let n = 0;
  for (let i = 0; i < total; i += stride) {
    const g = d[i * 4 + 1];
    if (g >= 254 || g <= 1) clipped++;
    n++;
  }
  return n ? clipped / n : 0;
}

/** Convenience: a canvas sized to the track's real frame, not the CSS box. */
export function makeFrameCanvas(track) {
  const { width, height } = track.getSettings();
  const canvas = document.createElement('canvas');
  canvas.width = width || 1280;
  canvas.height = height || 720;
  const ctx = canvas.getContext('2d', { willReadFrequently: true });
  return { canvas, ctx, width: canvas.width, height: canvas.height };
}
