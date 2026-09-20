# Four convincing measurements that were not measurements

Every number in this file was produced by this instrument, on a real MacBook
webcam, during its own development. Each one looked like a successful
measurement. None of them was.

We are publishing them because the research pass that preceded this build
concluded that a quantified account of what a browser camera path does to the
rolling-shutter signature does not exist in the measurement literature. We went
looking for a signal and found four ways to be fooled instead, which turns out
to be the more useful result.

The pattern worth taking away: **every gate we added was defeated by the next
scene we pointed at.** Confidence in this domain is not one number.

---

## 1. Peak over median, 394.7, on a blank wall

**What it showed:** a confidence ratio of 394.7 pointing at an ordinary wall
with no periodic source anywhere in frame.

**Why it happened:** the metric compared the tallest spectral bin to the
**global** median of the spectrum. A row profile of any real photograph is
dominated by low spatial frequency, because a bright ceiling and a dark desk
and lens vignetting are all slow gradients. Those low bins tower over the
high-bin noise floor, so the global ratio is enormous for any scene at all.

**What it means generally:** a confidence number that compares one part of a
spectrum to a different part of the same spectrum is measuring the spectrum's
shape, not the presence of a signal.

**Fix:** prominence against a running-median **local** baseline, so a bin is
compared to its own neighbourhood.

---

## 2. A four-row period, coherence 0.947, frame-to-frame spread 0.04%

**What it showed:** bin 512.63 of a 2048-point transform, which is a period of
almost exactly four rows. Cross-strip coherence 0.947. Frame-to-frame spread
0.04%.

Those are better numbers than any real measurement in this project has
produced. That is what made it dangerous.

**Why it happened:** no lamp has a four-row period. This was a fixed spatial
pattern from the imaging pipeline, most likely demosaic or rescaling structure.
A fixed pattern is **perfectly stable and perfectly coherent by construction**,
so the two gates designed to establish trustworthiness both scored it highest.

**What it means generally:** stability and spatial coherence are properties a
static artifact has *more* of than a real signal does. Tests that reward them
are inverted against this failure mode.

**Fix:** two. A ceiling in cycles per row, because above roughly 0.2 cycles per
row a "line" is a pattern repeating every few rows and that is the sensor, not
the world. And more fundamentally, background subtraction, below.

---

## 3. Bin 6.37, coherence 0.970, and physically consistent with 120 Hz

**What it showed:** 3.0555e-3 cycles per row, spread 1.31% over 14 frames,
motion 0.034, zero clipped pixels, pointed at a lit ceiling.

This one is the hardest, because it is **arithmetically plausible in every
respect**. If that reading were 120 Hz, the implied row time is 25.5 µs, which
would put 27.5 ms of readout inside a 33.3 ms frame with about 17% blanking.
That is exactly what a 1080p webcam should look like, and it sits between
horshack's dedicated-camera figures and the 34.38 µs that Sensors 26(16):5231
measured on an unmodified webcam.

**Why it happened:** the row profile showed one bump and one dip, then flat. A
3.3-cycle banding pattern must show three evenly spaced ripples across the full
height. An isolated excursion about one sixth of the frame tall puts its
spectral peak near bin 6 for reasons that have nothing to do with a lamp, so
bin position could not separate the two, and neither could coherence: a
horizontal band and a drifting horizontal edge both span every column.

**What it means generally:** a reading can pass a plausibility check against
independent published values and still be an artifact. Consistency with the
literature is not evidence of a measurement.

**Fix:** fit a sinusoid at the detected frequency by least squares and report
what fraction of the signal's energy it explains. A periodic line explains
nearly all of it; a localized bump explains very little, because its energy is
spread across many bins. Measured separation on the two shapes at the same
bin: **1.000 against 0.311.**

---

## 4. Coherence 0.995, caught only by the fit

**What it showed:** bin 4.08, cross-strip coherence **0.995**, motion 0.062,
0.6% clipped, pointed at a ceiling corner.

0.995 is as convincing as the coherence gate can possibly look.

**Why it was refused:** one wave at that frequency explained **39%** of the
trace. The instrument printed `not periodic` and drew the mismatch.

**What it means generally:** the gate that catches a given artifact is usually
not the gate you added last time. Four scenes, four different gates doing the
work.

---

## The gate that changed the most

Background subtraction, and it is a temporal argument rather than a spectral
one.

The scene, the vignetting and the sensor's own fixed pattern are identical in
every frame. Illumination flicker is not, because the frame period is not
locked to the mains, so its phase walks. A running mean of the row profile
therefore converges on the scene while the flicker averages itself away, and
the residual is the modulation.

Measured on a reconstruction of case 3 with the static artifact from case 2
deliberately made **stronger** than the lamp:

| | bin | prominence | coherence | outcome |
|---|---|---|---|---|
| raw frame | 6.79 | 1.77 | 0.703 | rejected |
| background subtracted | 6.63 | **34.21** | **1.000** | correct |

The four-row artifact does not get filtered out. It **disappears**, because it
is static.

---

## What we did not manage to measure

**We have no confirmed positive detection of mains flicker on this MacBook
webcam**, across a wall lamp, a lit ceiling, a ceiling corner and a desk scene.

We believe the cause is exposure, and it is structural rather than a defect in
this code. Rolling-shutter banding from a light only exists if a row's exposure
is shorter than one cycle of the source, which is 8.3 ms at 120 Hz. An indoor
webcam on auto-exposure typically sits at 16 to 33 ms, so the modulation is
integrated to zero before any code can see it.

The browser cannot correct this. `powerLineFrequency` is not a web constraint,
`exposureTime` exists only in the Image Capture spec as an optional advanced
constraint, iOS Safari does not implement Image Capture at all, and this
machine reports no exposure control of any kind.

That is why this instrument measures **rotation**, whose observable is
geometric and survives a long exposure, rather than illumination flicker, whose
observable is photometric and which the capture pipeline can erase before
delivery.

---

## Reproducing any of this

```
git clone https://github.com/StephenSook/rowclock
cd rowclock
node --test tests/*.test.mjs
```

No install, no dependencies, no network. Every case above is pinned as a
regression test so that none of them can be quietly tuned away later:

- case 1 as `REGRESSION: a scene with no periodic content must not score as confident`
- case 2 as `REGRESSION: a four-row pipeline artifact is excluded by the cycles-per-row ceiling`
- case 3 as `KEY DISCRIMINATOR: fitFraction is low for a single localized bump`
- the background-subtraction result as `BACKGROUND SUBTRACTION: a real line survives a stronger static artifact`
- and the single-frame limit itself as `DOCUMENTED LIMIT: a hard scene edge defeats the single-strip path at 120 Hz`
