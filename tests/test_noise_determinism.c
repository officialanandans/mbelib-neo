// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * Copyright (C) 2025 by arancormonk <180709949+arancormonk@users.noreply.github.com>
 */

/**
 * @file
 * @brief Determinism test for unvoiced noise synthesis using thread-local RNG.
 */
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include "mbelib-neo/mbelib.h"

/**
 * @brief Initialize deterministic unvoiced parameter sets for testing.
 * @param cur  Output current parameters (unvoiced bands populated).
 * @param prev Output previous parameters (copy of current).
 */
static void
fill_params_unvoiced(mbe_parms* cur, mbe_parms* prev) {
    mbe_parms prev_enh;                     // unused but required by init
    mbe_initMbeParms(cur, prev, &prev_enh); // resets and copies prev->cur

    // Overwrite cur/prev with deterministic unvoiced setup
    cur->w0 = 0.10f;
    cur->L = 24;
    for (int l = 1; l <= cur->L; ++l) {
        cur->Vl[l] = 0; // unvoiced
        cur->Ml[l] = 0.04f + 0.001f * (float)l;
        cur->PHIl[l] = 0.0f; // not used for unvoiced
        cur->PSIl[l] = 0.0f;
    }
    *prev = *cur; // start from same state
}

/**
 * @brief Test entry: verifies deterministic noise across identical seeds.
 */
int
main(void) {
    float out1[160], out2[160], out3[160];
    mbe_parms cur, prev;

    // First run with seed A
    mbe_setThreadRngSeed(0xA5A5A5u);
    fill_params_unvoiced(&cur, &prev);
    mbe_synthesizeSpeechf(out1, &cur, &prev, 8);

    // Second run with same seed A and same params -> identical buffer
    mbe_setThreadRngSeed(0xA5A5A5u);
    fill_params_unvoiced(&cur, &prev);
    mbe_synthesizeSpeechf(out2, &cur, &prev, 8);
    assert(memcmp(out1, out2, sizeof(out1)) == 0);

    // Third run with a different seed -> expect difference (very high probability)
    mbe_setThreadRngSeed(0x5A5A5Au);
    fill_params_unvoiced(&cur, &prev);
    mbe_synthesizeSpeechf(out3, &cur, &prev, 8);
    assert(memcmp(out1, out3, sizeof(out1)) != 0);
    return 0;
}
