"""Keyed random streams (DISCOVERY_LAYER_DESIGN §3.2): a draw is addressed by its key,
so appending candidates never moves the streams of existing ones. crc32 only turns a
string key into seed words; nothing is cached or reused by it."""

import zlib

import numpy as np


def stream(seed, *parts):
    words = [int(seed)] + [zlib.crc32(str(p).encode()) for p in parts]
    return np.random.default_rng(np.random.SeedSequence(words))
