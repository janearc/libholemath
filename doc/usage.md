# usage

    h := holemath.Hole{M: 1, A: 0.6}
    f := h.FrameOf(x, h.Static(x))      // a still observer's tetrad
    r := h.Arriving(x, f, n)            // light seen coming from n,
    end, fate := h.March(r, far, steps) //   marched back in time
    g := holemath.Shift(r.P, eye, end.P, emitter) // energy, eye over emitter
