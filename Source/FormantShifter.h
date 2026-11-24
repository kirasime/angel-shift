// Angel Shift - formant shifter core (no external dependencies).
//
// Idea: source-filter model with a causal, backward-adaptive LPC.
//   * The vocal-tract envelope is estimated from PAST samples only
//     -> zero algorithmic latency.
//   * The (pre-emphasised) input is whitened with A(z); the residual
//     (excitation) carries pitch and transients and is NOT touched
//     -> pitch is preserved, transients survive.
//   * The residual is re-synthesised with 1/A'(z), where A' is A with every
//     formant pole angle multiplied by the shift ratio (pole radii / bandwidths
//     stay) -> formants move, no new sharp resonances appear.
//   * The shift fades out above ~5 kHz, the top band is left alone.
//   * Noisy / unvoiced sounds, clicks and plosive onsets are detected with a
//     fast detector and stay dry (sample-accurate crossfade, still zero latency).
//   * At 0 semitones A'(z) == A(z), so the whole chain is an exact identity.
//   * Both filters are interpolated per sample between envelope updates.
#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>

namespace angel
{
class FormantShifter
{
public:
    static constexpr int kOrder = 32;   // LPC order (16 pole pairs)
    static constexpr int kHop   = 64;   // envelope update interval, full-rate samples

    void prepare (double sampleRate)
    {
        fs = sampleRate > 8000.0 ? sampleRate : 8000.0;

        decim = 1;   // full-rate analysis (a decimated model leaves the top band wildly amplified by A(z))
        fd    = fs / decim;

        lambdaD = std::exp (-1.0 / (0.012 * fd));            // 12 ms analysis window
        for (int k = 0; k <= kOrder; ++k)                     // lag window (~30 Hz)
        {
            const double x = 2.0 * kPi * 30.0 * k / fd;
            lagWin[k] = std::exp (-0.5 * x * x);
        }

        hp.designHP (std::min (4000.0, 0.45 * fs), fs, 0.7071);

        if (decim > 1)                                        // anti-alias filter, 4th-order Butterworth
        {
            lp[0].design (0.45 * fd, fs, 0.5412);
            lp[1].design (0.45 * fd, fs, 1.3066);
        }

        lambda     = std::exp (-1.0 / (0.012  * fs));
        lambdaLvl  = std::exp (-1.0 / (0.040  * fs));
        gAutoCoef  = 1.0 - std::exp (-1.0 / (0.015 * fs));
        peakDecay  = std::exp (-1.0 / (0.050  * fs));
        lambdaFast = std::exp (-1.0 / (0.0015 * fs));
        lambdaOnF  = std::exp (-1.0 / (0.00025 * fs));
        lambdaOnS  = std::exp (-1.0 / (0.008  * fs));
        onsetRel   = std::exp (-1.0 / (0.004  * fs));
        protectAtk = 1.0 - std::exp (-1.0 / (0.0003 * fs));
        protectRel = 1.0 - std::exp (-1.0 / (0.006  * fs));
        shiftCoef  = 1.0 - std::exp (-double (kHop) / (0.015 * fs));
        voicedAtk  = 1.0 - std::exp (-double (kHop) / (0.003 * fs));
        voicedRel  = 1.0 - std::exp (-double (kHop) / (0.020 * fs));
        thLo1      = 2.0 * kPi *  60.0 / fs;
        thLo2      = 2.0 * kPi * 150.0 / fs;

        reset();
    }

    void reset() noexcept
    {
        std::memset (r, 0, sizeof (r));
        ru0 = rf0 = ehfS = ehfF = efast = eslow = onset = 0.0;
        hp.clear();
        protect = 1.0;
        pin = pout = 0.0; gAuto = 1.0; pk[0] = pk[1] = 0.0;
        mPrev = 0.0;
        decCount = 0;
        for (auto& f : lp) f.clear();
        hs.clear();
        resetStates();

        setIdentity (aPrev); setIdentity (aNew);
        setIdentity (bPrev); setIdentity (bNew);
        gPrev = gNew = 1.0;
        byPrev = byNew = 1.0;
        hopPos = 0;
        warmValid = false;
        shiftSm = target;
        voiced = 1.0;
    }

    // Formant shift in semitones (-12 .. +12). Positive = formants up.
    void setShiftSemitones (float st) noexcept { target = (double) st; }

    // In-place processing. right may be nullptr (mono).
    void process (float* left, float* right, int numSamples) noexcept
    {
        const int nCh = right != nullptr ? 2 : 1;

        for (int n = 0; n < numSamples; ++n)
        {
            float* io[2] = { left + n, right != nullptr ? right + n : nullptr };

            const double in0 = (double) *io[0];
            const double in1 = nCh == 2 ? (double) *io[1] : in0;
            analyse (0.5 * (in0 + in1));

            // per-sample interpolation of the reflection coefficients
            // (|k| < 1 on both ends  =>  every intermediate filter is stable)
            const double t = (double) hopPos / (double) kHop;
            double ka[kOrder + 1], ks[kOrder + 1];
            for (int k = 1; k <= kOrder; ++k)
            {
                ka[k] = aPrev[k] + (aNew[k] - aPrev[k]) * t;
                ks[k] = bPrev[k] + (bNew[k] - bPrev[k]) * t;
            }
            const double g = gPrev + (gNew - gPrev) * t;
            const double byp = byPrev + (byNew - byPrev) * t;
            const double dry = std::max (std::max (protect, onset), byp);

            double eIn = 0.0, eOut = 0.0;
            for (int c = 0; c < nCh; ++c)
            {
                const double xr = (double) *io[c];
                const double x  = xr - kPre * ch[c].xRaw;   // pre-emphasis: removes the source tilt,
                ch[c].xRaw = xr;                            // so only the vocal-tract envelope is shifted

                // whitening A(z): lattice, forward error out, backward errors as state
                double* ba = ch[c].ba;
                double f = x, bcur = x;
                for (int m = 1; m <= kOrder; ++m)
                {
                    const double bOld = ba[m - 1];
                    const double fn = f + ka[m] * bOld;
                    const double bn = bOld + ka[m] * f;
                    ba[m - 1] = bcur;
                    bcur = bn;
                    f = fn;
                }

                // synthesis 1/A'(z): all-pole lattice
                double* bs = ch[c].bs;
                double bnew[kOrder + 1];
                double fc = f * g;                          // level compensation applied to the excitation
                for (int m = kOrder; m >= 1; --m)
                {
                    fc -= ks[m] * bs[m - 1];
                    bnew[m] = bs[m - 1] + ks[m] * fc;
                }
                double y = fc;
                bnew[0] = y;
                for (int m = 0; m < kOrder; ++m) bs[m] = bnew[m];

                if (! (std::abs (y) < 1.0e3 * (1.0 + std::abs (x))))   // NaN / blow-up guard
                {
                    resetStates();
                    y = x;
                }

                const double outRaw = y + kPre * ch[c].outPrev;    // de-emphasis (exact inverse at 0 st)
                ch[c].outPrev = outRaw;
                eIn  += xr * xr;
                eOut += outRaw * outRaw;

                // slow level match + envelope-following safety clamp (never relevant for normal material)
                double out = outRaw * gAuto;
                pk[c] = std::max (std::abs (xr), pk[c] * peakDecay);
                const double lim = 4.0 * pk[c] + 1.0e-4;
                out = std::clamp (out, -lim, lim);

                *io[c] = (float) (out + dry * (xr - out));             // noisy/transient parts stay dry
            }

            pin  = lambdaLvl * pin  + (1.0 - lambdaLvl) * eIn  / nCh;
            pout = lambdaLvl * pout + (1.0 - lambdaLvl) * eOut / nCh;
            const double gT = std::clamp (std::sqrt ((pin + 1.0e-12) / (pout + 1.0e-12)), 0.1, 10.0);
            gAuto += (gT - gAuto) * gAutoCoef;

            if (++hopPos == kHop)
            {
                hopPos = 0;
                std::memcpy (aPrev, aNew, sizeof (aPrev));
                std::memcpy (bPrev, bNew, sizeof (bPrev));
                gPrev = gNew;
                byPrev = byNew;
                updateFilters();
            }
        }
    }

    static constexpr int latencySamples() noexcept { return 0; }

private:
    static constexpr double kPi   = 3.14159265358979323846;
    static constexpr double kRMax = 0.992;      // max pole radius = min. formant bandwidth ~ 120 Hz @48k (as narrow as a natural vowel, no ringing)
    static constexpr double kPre  = 0.9;        // pre/de-emphasis coefficient

    using Cx = std::complex<double>;

    struct Hist
    {
        double buf[2 * kOrder] = {};
        int pos = 0;
        void clear() noexcept { std::memset (buf, 0, sizeof (buf)); pos = 0; }
        void push (double v) noexcept
        {
            pos = pos > 0 ? pos - 1 : kOrder - 1;
            buf[pos] = buf[pos + kOrder] = v;
        }
        const double* last() const noexcept { return buf + pos; }   // last()[0] = newest
    };
    struct Channel { double ba[kOrder] = {}, bs[kOrder] = {}, xRaw = 0.0, outPrev = 0.0; };

    struct Biquad
    {
        double b0 = 1, b1 = 0, b2 = 0, a1 = 0, a2 = 0, z1 = 0, z2 = 0;
        void clear() noexcept { z1 = z2 = 0.0; }
        void design (double fc, double sr, double q) noexcept
        {
            const double w = 2.0 * kPi * fc / sr, cs = std::cos (w), al = std::sin (w) / (2.0 * q);
            const double a0 = 1.0 + al;
            b0 = (1.0 - cs) * 0.5 / a0; b1 = (1.0 - cs) / a0; b2 = b0;
            a1 = -2.0 * cs / a0;        a2 = (1.0 - al) / a0;
        }
        void designHP (double fc, double sr, double q) noexcept
        {
            const double w = 2.0 * kPi * fc / sr, cs = std::cos (w), al = std::sin (w) / (2.0 * q);
            const double a0 = 1.0 + al;
            b0 = (1.0 + cs) * 0.5 / a0; b1 = -(1.0 + cs) / a0; b2 = b0;
            a1 = -2.0 * cs / a0;         a2 = (1.0 - al) / a0;
        }
        double process (double x) noexcept
        {
            const double y = b0 * x + z1;
            z1 = b1 * x - a1 * y + z2;
            z2 = b2 * x - a2 * y;
            return y;
        }
    };

    static void setIdentity (double* a) noexcept
    {
        std::memset (a, 0, sizeof (double) * (kOrder + 1));
        a[0] = 1.0;
    }

    void resetStates() noexcept
    {
        for (auto& c : ch) { std::memset (c.ba, 0, sizeof (c.ba)); std::memset (c.bs, 0, sizeof (c.bs)); c.xRaw = 0.0; c.outPrev = 0.0; }
    }

    static double smoothstep (double e0, double e1, double v) noexcept
    {
        const double t = std::clamp ((v - e0) / (e1 - e0), 0.0, 1.0);
        return t * t * (3.0 - 2.0 * t);
    }

    // ---- analysis (per sample, O(order)) ----------------------------------
    void analyse (double m) noexcept
    {
        // LPC input: pre-emphasised, band-limited and decimated copy of the mid signal
        double s = m - kPre * mPrev;
        if (decim > 1) { s = lp[0].process (s); s = lp[1].process (s); }
        if (++decCount >= decim)
        {
            decCount = 0;
            const double* sp = hs.last();
            r[0] = lambdaD * r[0] + s * s;
            for (int k = 1; k <= kOrder; ++k) r[k] = lambdaD * r[k] + s * sp[k - 1];
            hs.push (s);
        }

        // noisiness cue: fraction of the energy above 4 kHz (fricatives ~ 0.5+, vowels / whispered vowels ~ < 0.15)
        const double h = hp.process (m);
        ru0  = lambda     * ru0  + m * m;
        ehfS = lambda     * ehfS + h * h;          // 12 ms
        rf0  = lambdaFast * rf0  + m * m;
        ehfF = lambdaFast * ehfF + h * h;          // 1.5 ms: sample-accurate protection

        const double hfF = ehfF / (rf0 + 1.0e-12);
        const double pT = smoothstep (0.25, 0.6, hfF);
        protect += (pT - protect) * (pT > protect ? protectAtk : protectRel);

        // onset detector: sudden jump of the first-difference energy (clicks, plosives)
        const double d = m - mPrev;
        efast = lambdaOnF * efast + (1.0 - lambdaOnF) * d * d;
        eslow = lambdaOnS * eslow + (1.0 - lambdaOnS) * d * d;
        const double ratio = efast / (eslow + 1.0e-12);
        const double oT = efast > 1.0e-8 ? smoothstep (8.0, 20.0, ratio) : 0.0;
        onset = std::max (onset * onsetRel, oT);

        mPrev = m;
    }

    // ---- Levinson-Durbin ---------------------------------------------------
    static bool levinson (const double* rr, double* a) noexcept
    {
        double tmp[kOrder + 1];
        setIdentity (a);
        double err = rr[0];
        if (! (err > 0.0)) return false;
        for (int i = 1; i <= kOrder; ++i)
        {
            double acc = rr[i];
            for (int j = 1; j < i; ++j) acc += a[j] * rr[i - j];
            const double k = -acc / err;
            if (! (std::abs (k) < 0.99999)) return false;
            std::memcpy (tmp, a, sizeof (tmp));
            for (int j = 1; j < i; ++j) a[j] = tmp[j] + k * tmp[i - j];
            a[i] = k;
            err *= (1.0 - k * k);
            if (! (err > 0.0)) return false;
        }
        return true;
    }

    // power gain of 1/A(z) for unit white noise: prod 1/(1-k^2); -1 if unstable
    static double arPowerGain (const double* a, int order, double* kOut = nullptr) noexcept
    {
        double cur[kOrder + 2], nxt[kOrder + 2];
        std::memcpy (cur, a, sizeof (double) * (order + 1));
        double g = 1.0;
        for (int m = order; m >= 1; --m)
        {
            const double k = cur[m];
            const double d = 1.0 - k * k;
            if (! (std::abs (k) < 0.9999) || ! (d > 1.0e-12)) return -1.0;
            g /= d;
            if (kOut != nullptr) kOut[m] = k;
            for (int i = 1; i < m; ++i) nxt[i] = (cur[i] - k * cur[m - i]) / d;
            for (int i = 1; i < m; ++i) cur[i] = nxt[i];
        }
        return g;
    }

    // multiply poly (a[0]=1, degree deg) by (1 + c1 z^-1 [+ c2 z^-2]); returns new degree
    static int mulFactor (double* a, int deg, double c1, double c2, bool quad) noexcept
    {
        const int nd = deg + (quad ? 2 : 1);
        for (int i = nd; i >= 1; --i)
        {
            double v = a[i] + c1 * a[i - 1];
            if (quad && i >= 2) v += c2 * a[i - 2];
            a[i] = v;
        }
        return nd;
    }

    // power gain of the complete synthesis chain 1/(A(z) (1 - pre z^-1)), evaluated at the decimated rate
    // (the model only describes 0-8 kHz; computing it at the full rate would be dominated by the empty top band)
    static double chainPowerGain (const double* a, double pre) noexcept
    {
        double t[kOrder + 2] = {};
        std::memcpy (t, a, sizeof (double) * (kOrder + 1));
        mulFactor (t, kOrder, -pre, 0.0, false);
        return arPowerGain (t, kOrder + 1);
    }

    // ---- root finding (Aberth-Ehrlich, warm-started) ----------------------
    bool findRoots (const double* a) noexcept
    {
        constexpr int p = kOrder;
        if (! warmValid)
            for (int i = 0; i < p; ++i)
                rt[i] = std::polar (0.7 + 0.25 * double (i) / p, 2.0 * kPi * (i + 0.5) / p + 0.4);

        bool ok = false;
        for (int it = 0; it < 80; ++it)
        {
            double worst = 0.0;
            for (int i = 0; i < p; ++i)
            {
                const Cx z = rt[i];
                Cx pv (a[0], 0.0), dv (0.0, 0.0);
                for (int k = 1; k <= p; ++k) { dv = dv * z + pv; pv = pv * z + a[k]; }
                if (std::abs (dv) < 1.0e-300) { rt[i] = z + Cx (1.0e-6, 1.0e-6); worst = 1.0; continue; }
                const Cx w = pv / dv;
                Cx sum (0.0, 0.0);
                for (int j = 0; j < p; ++j)
                    if (j != i) sum += 1.0 / (z - rt[j]);
                const Cx den = 1.0 - w * sum;
                const Cx step = std::abs (den) > 1.0e-300 ? w / den : w;
                rt[i] = z - step;
                worst = std::max (worst, std::abs (step));
            }
            if (worst < 1.0e-9) { ok = true; break; }
            if (it == 79 && worst < 1.0e-6) ok = true;
        }

        if (ok)
            for (int i = 0; i < p; ++i)
                if (! std::isfinite (rt[i].real()) || ! std::isfinite (rt[i].imag())) ok = false;

        warmValid = ok;
        return ok;
    }

    // ---- envelope update (every kHop samples) -----------------------------
    void updateFilters() noexcept
    {
        // parameter smoothing
        shiftSm += (target - shiftSm) * shiftCoef;

        const double hfS = ehfS / (ru0 + 1.0e-12);
        const double vt = 1.0 - smoothstep (0.25, 0.6, hfS);
        voiced += (vt - voiced) * (vt < voiced ? voicedAtk : voicedRel);
        const double st = shiftSm * voiced;           // noisy sounds stay unshifted
        byNew = 1.0 - smoothstep (0.02, 0.3, std::abs (st));   // centred knob = exact bypass

        if (! (r[0] > 1.0e-10)) return;                // silence: hold the filters

        double rr[kOrder + 1];
        for (int k = 0; k <= kOrder; ++k) rr[k] = r[k] * lagWin[k];
        rr[0] = rr[0] * (1.0 + 3.0e-4) + 1.0e-12;

        double a[kOrder + 1];
        if (! levinson (rr, a)) return;
        if (! findRoots (a)) return;

        // map the poles from the decimated rate back to the full rate
        double pr[kOrder / 2], pth[kOrder / 2], reals[kOrder];
        int np = 0, nr = 0;
        const double invD = 1.0 / decim;
        for (int i = 0; i < kOrder; ++i)
        {
            const double re = rt[i].real(), im = rt[i].imag();
            if (std::abs (im) < 1.0e-7)
            {
                if (re > 0.0)                                  // negative real poles (Nyquist_dec) are dropped
                    reals[nr++] = std::min (std::pow (re, invD), kRMax);
            }
            else if (im > 0.0 && np < kOrder / 2)
            {
                pr[np]  = std::min (std::pow (std::hypot (re, im), invD), kRMax);
                pth[np] = std::atan2 (im, re) * invD;
                ++np;
            }
        }

        // A_c: envelope with radii clamped (used for whitening AND as the st=0 synthesis)
        double ac[kOrder + 2] = {}, acD[kOrder + 2] = {};
        ac[0] = acD[0] = 1.0;
        int deg = 0, degAD = 0;
        for (int i = 0; i < np; ++i)
        {
            deg   = mulFactor (ac,  deg,   -2.0 * pr[i] * std::cos (pth[i]), pr[i] * pr[i], true);
            degAD = mulFactor (acD, degAD, -2.0 * std::pow (pr[i], decim) * std::cos (pth[i] * decim), std::pow (pr[i], 2 * decim), true);
        }
        for (int i = 0; i < nr; ++i)
        {
            deg   = mulFactor (ac,  deg,   -reals[i], 0.0, false);
            degAD = mulFactor (acD, degAD, -std::pow (reals[i], decim), 0.0, false);
        }

        // A': formant angles scaled, radii kept (soft fades at the band edges)
        const double ratio = std::exp2 (st / 12.0);
        const double c     = std::clamp (std::abs (st) / 0.5, 0.0, 1.0);   // fades engage smoothly from 0
        double bc[kOrder + 2] = {}, bcD[kOrder + 2] = {};
        bc[0] = bcD[0] = 1.0;
        deg = 0;
        int degBD = 0;
        for (int i = 0; i < np; ++i)
        {
            // Formants live below ~5 kHz; the shift is faded out towards the top of the model band.
            const double hz = pth[i] * fs / (2.0 * kPi);
            const double w  = 1.0 - smoothstep (4500.0, 8000.0, hz);
            double th = pth[i] * (1.0 + (ratio - 1.0) * w);
            const double f1 = std::clamp ((0.98 * kPi - th) / (0.08 * kPi), 0.0, 1.0); // fade near Nyquist
            const double f2 = smoothstep (thLo1, thLo2, th);                            // fade near DC (no rumble)
            th = std::min (th, 0.999 * kPi);
            const double rad = pr[i] * (1.0 + c * (f1 * f2 - 1.0));
            deg = mulFactor (bc, deg, -2.0 * rad * std::cos (th), rad * rad, true);
            if (th * decim < 0.98 * kPi)                    // poles pushed beyond the modelled band leave the gain budget
                degBD = mulFactor (bcD, degBD, -2.0 * std::pow (rad, decim) * std::cos (th * decim), std::pow (rad, 2 * decim), true);
        }
        for (int i = 0; i < nr; ++i)
        {
            deg   = mulFactor (bc,  deg,   -reals[i], 0.0, false);
            degBD = mulFactor (bcD, degBD, -std::pow (reals[i], decim), 0.0, false);
        }

        // reflection coefficients of both filters (stable by construction: radii <= kRMax)
        double kAc[kOrder + 2] = {}, kBc[kOrder + 2] = {};
        if (! (arPowerGain (ac, kOrder, kAc) > 0.0)) return;
        if (! (arPowerGain (bc, kOrder, kBc) > 0.0))      // unstable result: fall back to "no shift"
        {
            std::memcpy (bc, ac, sizeof (bc));
            std::memcpy (kBc, kAc, sizeof (kBc));
        }

        // loudness compensation from AR power gains (including the de-emphasis)
        double g = 1.0;
        const double preD = std::pow (kPre, decim);
        const double po = chainPowerGain (acD, preD), pn = chainPowerGain (bcD, preD);
        if (po > 0.0 && pn > 0.0) g = std::clamp (std::sqrt (po / pn), 0.005, 200.0);   // first guess; the slow level follower trims the rest

        std::memcpy (aNew, kAc, sizeof (double) * (kOrder + 1));
        std::memcpy (bNew, kBc, sizeof (double) * (kOrder + 1));
        gNew = g;
    }

    // ---- state -------------------------------------------------------------
    double fs = 48000.0, fd = 16000.0;
    int decim = 3, decCount = 0;
    double lambda = 0.99, lambdaD = 0.99, lambdaFast = 0.9, lambdaOnF = 0.8, lambdaOnS = 0.99, onsetRel = 0.99;
    double shiftCoef = 0.1, voicedAtk = 0.5, voicedRel = 0.1, protectAtk = 0.1, protectRel = 0.01;
    double thLo1 = 0.0, thLo2 = 0.0;
    double lagWin[kOrder + 1] = {};
    Biquad lp[2];

    double r[kOrder + 1] = {};
    double ru0 = 0.0, rf0 = 0.0, ehfS = 0.0, ehfF = 0.0, mPrev = 0.0;
    Biquad hp;
    double efast = 0.0, eslow = 0.0, onset = 0.0, protect = 1.0;
    double pin = 0.0, pout = 0.0, gAuto = 1.0, lambdaLvl = 0.99, gAutoCoef = 0.01, peakDecay = 0.99;
    double pk[2] = { 0.0, 0.0 };
    Hist hs;
    Channel ch[2];

    double aPrev[kOrder + 1] = {}, aNew[kOrder + 1] = {};
    double bPrev[kOrder + 1] = {}, bNew[kOrder + 1] = {};
    double gPrev = 1.0, gNew = 1.0, byPrev = 1.0, byNew = 1.0;
    int hopPos = 0;

    Cx rt[kOrder];
    bool warmValid = false;

    double target = 0.0, shiftSm = 0.0, voiced = 1.0;
};
} // namespace angel
